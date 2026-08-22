#!/bin/bash

################################################################################
#                                                                              #
#               GridCAT Analysis Pipeline — Interactive Master Runner          #
#                                                                              #
#  WHAT THIS SCRIPT DOES:                                                     #
#  ──────────────────────                                                     #
#  This is the single entry point for running the complete GridCAT            #
#  neuroimaging analysis pipeline. It provides an interactive menu that       #
#  guides users through each analysis step, handling SLURM job submission,     #
#  configuration validation, and results collection.                          #
#                                                                              #
#  HOW TO USE:                                                                #
#  ──────────                                                                 #
#  1. Edit pipeline_config.cfg in the same directory as this script to        #
#     set your paths, parameters, and HPC settings                            #
#  2. Run: bash run_pipeline.sh                                               #
#  3. Follow the interactive menu to run steps                                #
#                                                                              #
#  WHAT IS pipeline_config.cfg:                                              #
#  ───────────────────────────                                                #
#  The configuration file (pipeline_config.cfg) contains all settings for     #
#  the pipeline: data paths, preprocessing parameters, HPC cluster settings,  #
#  and study design parameters. You MUST edit this file before running the    #
#  pipeline. The master script will source this file automatically.           #
#                                                                              #
#  PIPELINE STEPS:                                                            #
#  ───────────────                                                            #
#  Step 0: Discover subjects & sessions from BIDS directory                   #
#          Scans the BIDS folder and creates a list of subjects/sessions      #
#          to analyze. (Runs step0_make_subses_list.sh)                       #
#                                                                              #
#  Step 1: Validate configuration                                             #
#          Checks that all required files, directories, and settings exist    #
#          and are valid. (Runs validate_config.sh)                           #
#                                                                              #
#  Step 2: Run SPM preprocessing (SLURM)                                      #
#          Submits preprocessing jobs to the HPC cluster as a SLURM array.    #
#          Each subject/session is processed in parallel.                      #
#          (Submits spm_preproc_array.sbatch)                                 #
#                                                                              #
#          Step 2 is itself split into five stages that can be run on their   #
#          own: topup -> vdm -> realign -> coreg -> smooth. Pick them with     #
#          menu option 2S, set PREPROC_STAGES in the config, or run            #
#          submit_preproc.sh --stages <list> outside this menu. Use this to    #
#          resume from fieldmaps you produced yourself, or to redo just the    #
#          smoothing without repeating the realignment.                        #
#                                                                              #
#  Step 3: Prepare data for GridCAT                                           #
#          Reorganizes preprocessed data into the format required by GridCAT. #
#          (Runs run_prep.sh which calls prepare_gridcat_directory.m)         #
#                                                                              #
#  Step 4: Run GridCAT analysis (SLURM)                                       #
#          Submits the GridCAT GLM statistical analysis to HPC.               #
#          (Submits run_gridcat.sh which calls run_gridcat_analysis.m)        #
#                                                                              #
#  TYPICAL WORKFLOW:                                                           #
#  ────────────────                                                            #
#  • Run "ALL steps" (option A) once. It runs the local steps (0,              #
#    filter, 1) right away, then submits the three cluster stages              #
#    (2 → 3 → 4) as a single SLURM dependency chain.                           #
#  • SLURM holds each stage until the previous one finishes, so the            #
#    whole pipeline runs unattended — no relaunching between stages.           #
#  • Use "Check job status" (option C) to watch the chain progress.            #
#  • Run a single stage on its own by picking its number.                      #
#                                                                              #
#  TIP FOR NEW USERS:                                                         #
#  ──────────────────                                                         #
#  If you're unsure, start with "Show current settings" (option S) to verify  #
#  your configuration is correct. Use "Check job status" (option C) to see    #
#  if your submitted jobs are running on the cluster.                         #
#                                                                              #
################################################################################

set -euo pipefail

# Color codes for formatted output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'  # No Color

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.cfg"
LOG_DIR="${SCRIPT_DIR}/logs"

# --- Run-all chaining state ---
# ASSUME_YES      : when "true", ask_confirmation auto-approves (used by run-all).
# NEXT_DEPENDENCY : SLURM dependency spec applied to the next submitted job
#                   (e.g. "afterok:1234"); empty for a standalone submission.
# LAST_JOB_ID     : job id of the most recent successful submission (set by steps
#                   2/3/4 so run-all can chain the next stage onto it).
# PREPROC_STAGES_OVERRIDE : stage list chosen interactively for one Step 2
#                   submission; empty means use PREPROC_STAGES from the config.
ASSUME_YES=false
NEXT_DEPENDENCY=""
LAST_JOB_ID=""
PREPROC_STAGES_OVERRIDE=""

# Create logs directory if it doesn't exist
mkdir -p "${LOG_DIR}"

# Initialize log file with timestamp
LOG_FILE="${LOG_DIR}/pipeline_run_$(date +%Y%m%d_%H%M%S).log"

################################################################################
# UTILITY FUNCTIONS
################################################################################

# Log a message to both stdout and log file
log_message() {
    local level=$1
    shift
    local message="$@"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')

    echo "[${timestamp}] [${level}] ${message}" >> "${LOG_FILE}"
    echo -e "${message}"
}

# Print colored header
print_header() {
    local title="$1"
    printf "\n"
    printf "${BLUE}╔════════════════════════════════════════════════════╗${NC}\n"
    printf "${BLUE}║${NC} %-48s ${BLUE}║${NC}\n" "$title"
    printf "${BLUE}╚════════════════════════════════════════════════════╝${NC}\n"
    printf "\n"
}

# Print colored success message
print_success() {
    printf "${GREEN}✓ SUCCESS:${NC} $@\n"
}

# Print colored error message
print_error() {
    printf "${RED}✗ ERROR:${NC} $@\n"
}

# Print colored warning message
print_warning() {
    printf "${YELLOW}⚠ WARNING:${NC} $@\n"
}

# Print colored info message
print_info() {
    printf "${BLUE}ℹ INFO:${NC} $@\n"
}

# Ask user for confirmation (Y/N)
# In run-all mode (ASSUME_YES=true) this auto-approves without prompting.
ask_confirmation() {
    local prompt="$1"
    local response

    if [[ "${ASSUME_YES}" == "true" ]]; then
        printf "${prompt} (yes/no): yes  ${BLUE}[auto]${NC}\n"
        return 0
    fi

    while true; do
        printf "${prompt} (yes/no): "
        read -r response
        case "${response}" in
            [Yy][Ee][Ss]|[Yy])
                return 0
                ;;
            [Nn][Oo]|[Nn])
                return 1
                ;;
            *)
                printf "Please answer 'yes' or 'no'.\n"
                ;;
        esac
    done
}

# Check if config file exists
check_config() {
    if [[ ! -f "${CONFIG_FILE}" ]]; then
        print_error "Configuration file not found: ${CONFIG_FILE}"
        printf "\nPlease create the file '${CONFIG_FILE}' with your pipeline settings.\n"
        return 1
    fi
    return 0
}

# Source the configuration file
source_config() {
    if ! source "${CONFIG_FILE}" 2>/dev/null; then
        print_error "Failed to source configuration file"
        return 1
    fi
    return 0
}

################################################################################
# STEP FUNCTIONS
################################################################################

# Step 0: Discover subjects and sessions
run_step_0() {
    print_header "STEP 0: Discover Subjects & Sessions"

    if [[ ! -f "${SCRIPT_DIR}/step0_make_subses_list.sh" ]]; then
        print_error "step0_make_subses_list.sh not found in ${SCRIPT_DIR}"
        return 1
    fi

    print_info "This step will scan the BIDS directory and create a list of subjects/sessions."
    printf "  Input:  ${BIDS_ROOT}\n"
    printf "  Output: ${SCRIPT_DIR}/subses_list.txt\n\n"

    if ! ask_confirmation "Run Step 0?"; then
        print_info "Step 0 skipped."
        return 0
    fi

    printf "\n"
    if bash "${SCRIPT_DIR}/step0_make_subses_list.sh"; then
        print_success "Step 0 completed successfully."
        log_message "INFO" "Step 0 completed successfully"
        return 0
    else
        print_error "Step 0 failed. Check the error message above."
        log_message "ERROR" "Step 0 failed"
        return 1
    fi
}

# Step 0b: Filter subject-session list
run_step_filter() {
    print_header "FILTER: Apply Subject/Session Filters"

    if [[ ! -f "${SCRIPT_DIR}/filter_subses_list.sh" ]]; then
        print_error "filter_subses_list.sh not found in ${SCRIPT_DIR}"
        return 1
    fi

    local SUBSES_LIST="${SCRIPT_DIR}/subses_list.txt"
    if [[ ! -f "${SUBSES_LIST}" ]]; then
        print_error "Subject/session list not found: ${SUBSES_LIST}"
        print_info "Please run Step 0 first to generate this file."
        return 1
    fi

    # Check if any filters are configured
    local has_filter=false
    if [[ -n "${INCLUDE_SUBJECTS:-}" ]] || [[ -n "${INCLUDE_SESSIONS:-}" ]] || \
       [[ -n "${EXCLUDE_SUBJECTS:-}" ]] || [[ -n "${EXCLUDE_SESSIONS:-}" ]] || \
       [[ "${REQUIRE_ROI:-false}" == "true" ]]; then
        has_filter=true
    fi

    if [[ "$has_filter" == "false" ]]; then
        print_info "No filters are configured in pipeline_config.cfg."
        printf "  INCLUDE_SUBJECTS: (empty)\n"
        printf "  INCLUDE_SESSIONS: (empty)\n"
        printf "  EXCLUDE_SUBJECTS: (empty)\n"
        printf "  EXCLUDE_SESSIONS: (empty)\n"
        printf "  REQUIRE_ROI:      false\n\n"
        print_info "Skipping filter — all subject-sessions will be processed."
        return 0
    fi

    print_info "Applying filters from pipeline_config.cfg:"
    [[ -n "${INCLUDE_SUBJECTS:-}" ]] && printf "  INCLUDE_SUBJECTS: ${INCLUDE_SUBJECTS}\n"
    [[ -n "${INCLUDE_SESSIONS:-}" ]] && printf "  INCLUDE_SESSIONS: ${INCLUDE_SESSIONS}\n"
    [[ -n "${EXCLUDE_SUBJECTS:-}" ]] && printf "  EXCLUDE_SUBJECTS: ${EXCLUDE_SUBJECTS}\n"
    [[ -n "${EXCLUDE_SESSIONS:-}" ]] && printf "  EXCLUDE_SESSIONS: ${EXCLUDE_SESSIONS}\n"
    [[ "${REQUIRE_ROI:-false}" == "true" ]] && printf "  REQUIRE_ROI:      true\n"
    printf "\n"

    if ! ask_confirmation "Apply these filters?"; then
        print_info "Filter step skipped."
        return 0
    fi

    printf "\n"
    if bash "${SCRIPT_DIR}/filter_subses_list.sh"; then
        print_success "Filter completed successfully."
        log_message "INFO" "Subject-session filter completed successfully"
        return 0
    else
        print_error "Filter failed. Check the error message above."
        log_message "ERROR" "Subject-session filter failed"
        return 1
    fi
}

# Step 1: Validate configuration
run_step_1() {
    print_header "STEP 1: Validate Configuration"

    if [[ ! -f "${SCRIPT_DIR}/validate_config.sh" ]]; then
        print_error "validate_config.sh not found in ${SCRIPT_DIR}"
        return 1
    fi

    print_info "This step will validate your configuration and check for required files/directories."
    printf "\n"

    if ! ask_confirmation "Run Step 1?"; then
        print_info "Step 1 skipped."
        return 0
    fi

    printf "\n"
    if bash "${SCRIPT_DIR}/validate_config.sh"; then
        print_success "Step 1 completed successfully."
        log_message "INFO" "Step 1 completed successfully"
        return 0
    else
        print_error "Step 1 failed. Check the error messages above."
        log_message "ERROR" "Step 1 failed"
        return 1
    fi
}

# Step 2: Run SPM preprocessing
run_step_2() {
    print_header "STEP 2: Run SPM Preprocessing (SLURM)"

    if [[ ! -f "${SCRIPT_DIR}/spm_preproc_array.sbatch" ]]; then
        print_error "spm_preproc_array.sbatch not found in ${SCRIPT_DIR}"
        return 1
    fi

    local SUBSES_LIST="${SCRIPT_DIR}/subses_list.txt"
    if [[ ! -f "${SUBSES_LIST}" ]]; then
        print_error "Subject/session list not found: ${SUBSES_LIST}"
        print_info "Please run Step 0 first to generate this file."
        return 1
    fi

    # Count the number of subjects/sessions
    local num_subjects
    num_subjects=$(wc -l < "${SUBSES_LIST}")

    # Stage list: an interactive choice wins over PREPROC_STAGES in the config.
    local stages="${PREPROC_STAGES_OVERRIDE:-${PREPROC_STAGES:-all}}"
    stages="${stages// /}"

    # Partial runs get a job name that says which stages they cover, so they
    # are recognisable in squeue next to a full preprocessing array.
    local job_name="spm_preproc"
    [[ "${stages}" != "all" ]] && job_name="spm_${stages//,/_}"

    print_info "This step will submit preprocessing jobs to the HPC cluster."
    printf "  Script:    spm_preproc_array.sbatch\n"
    printf "  Stages:    ${stages}\n"
    printf "  Mode:      ${PREPROC_MODE:-<not set>}\n"
    printf "  Array size: 1-%d (one job per subject/session)\n" "${num_subjects}"
    printf "  Subjects/sessions: %d\n" "${num_subjects}"
    printf "\nThese jobs will run in parallel on the cluster.\n"
    printf "Processing time depends on your data size and cluster availability.\n\n"

    if ! ask_confirmation "Submit Step 2 preprocessing jobs?"; then
        print_info "Step 2 skipped."
        return 0
    fi

    printf "\n"
    print_info "Submitting SLURM array job..."

    local -a dep_flag=()
    [[ -n "${NEXT_DEPENDENCY}" ]] && dep_flag=(--dependency="${NEXT_DEPENDENCY}")

    # Email once when the WHOLE array finishes. Omitting ARRAY_TASKS from
    # --mail-type makes SLURM treat the array as a single unit, so END/FAIL
    # generate one message for all subjects rather than one per task.
    local -a mail_flag=()
    if [[ -n "${NOTIFY_EMAIL:-}" ]]; then
        mail_flag=(--mail-user="${NOTIFY_EMAIL}" --mail-type=END,FAIL)
    fi
    LAST_JOB_ID=""

    local sbatch_output
    if sbatch_output=$(sbatch \
        --export=ALL,PIPELINE_SCRIPT_DIR="${SCRIPT_DIR}" \
        --job-name="${job_name}" \
        --partition="${SLURM_PARTITION}" \
        --time="${PREPROC_TIME}" \
        --cpus-per-task="${PREPROC_CPUS}" \
        --mem-per-cpu="${PREPROC_MEM_PER_CPU}" \
        --array=1-"${num_subjects}" \
        "${dep_flag[@]}" \
        "${mail_flag[@]}" \
        "${SCRIPT_DIR}/spm_preproc_array.sbatch" --stages "${stages}" 2>&1); then
        print_success "Preprocessing jobs submitted!"

        # Extract job ID from sbatch output
        local job_id
        job_id=$(echo "${sbatch_output}" | grep -oP 'Submitted batch job \K[0-9]+' || true)
        LAST_JOB_ID="${job_id}"

        if [[ -n "${job_id}" ]]; then
            printf "\n${GREEN}Job ID: ${job_id}${NC}\n"
            printf "Array Range: 1-%d\n" "${num_subjects}"
            printf "\nYou can check the status of your jobs using:\n"
            printf "  squeue -j ${job_id}\n\n"
            log_message "INFO" "Step 2 preprocessing jobs submitted with Job ID: ${job_id} (stages: ${stages})"
        else
            printf "\n${GREEN}${sbatch_output}${NC}\n"
            log_message "INFO" "Step 2 preprocessing jobs submitted"
        fi

        printf "Please wait for these jobs to complete on the cluster before running Step 3.\n"
        printf "You can check progress using the 'Check job status' menu option.\n"
        return 0
    else
        print_error "Failed to submit SLURM jobs"
        printf "Error output:\n${sbatch_output}\n"
        log_message "ERROR" "Step 2 failed to submit jobs: ${sbatch_output}"
        return 1
    fi
}

# Step 2 (partial): pick which preprocessing stages to submit
#
# The SPM preprocessing is split into five stages that each read what they need
# from disk, so any suffix of the pipeline can be re-run on its own. This lets
# you stop after the fieldmaps, swap in fieldmaps you made yourself, or redo
# just the smoothing without repeating the expensive realignment.
run_step_2_stages() {
    print_header "STEP 2 (partial): Choose Preprocessing Stages"

    printf "The preprocessing runs as five stages, in this order:\n\n"
    printf "  ${BLUE}topup${NC}    Estimate the distortion field with FSL topup\n"
    printf "           ${BLUE}->${NC} func/topup_results_*, func/topup_acqparams.txt\n"
    printf "  ${BLUE}vdm${NC}      Build one voxel displacement map per task\n"
    printf "           ${BLUE}->${NC} func/vdm_task-<N>.nii\n"
    printf "  ${BLUE}realign${NC}  Apply the correction + motion correction\n"
    printf "           ${BLUE}->${NC} func/u*_bold.nii, rp_*.txt, meanu_session.nii\n"
    printf "  ${BLUE}coreg${NC}    Coregister T2w -> session mean, reslice ROI masks\n"
    printf "           ${BLUE}->${NC} anat/r*_mask.nii\n"
    printf "  ${BLUE}smooth${NC}   Gaussian smoothing of the preprocessed BOLD\n"
    printf "           ${BLUE}->${NC} func/su*_bold.nii\n\n"

    printf "Presets:\n"
    printf "  ${BLUE}[1]${NC} all                        — the complete preprocessing\n"
    printf "  ${BLUE}[2]${NC} fieldmap                   — topup + vdm, then stop\n"
    printf "  ${BLUE}[3]${NC} vdm,realign,coreg,smooth   — resume from topup output made elsewhere\n"
    printf "  ${BLUE}[4]${NC} realign,coreg,smooth       — resume from ready-made VDMs\n"
    printf "  ${BLUE}[5]${NC} coreg,smooth               — redo coregistration and smoothing\n"
    printf "  ${BLUE}[6]${NC} smooth                     — re-smooth only\n"
    printf "  ${BLUE}[C]${NC} custom                     — type your own comma-separated list\n"
    printf "  ${BLUE}[X]${NC} cancel\n\n"

    printf "Enter your choice: "
    local pick
    read -r pick

    local stages=""
    case "${pick}" in
        1) stages="all" ;;
        2) stages="fieldmap" ;;
        3) stages="vdm,realign,coreg,smooth" ;;
        4) stages="realign,coreg,smooth" ;;
        5) stages="coreg,smooth" ;;
        6) stages="smooth" ;;
        [Cc])
            printf "Stages (comma-separated, e.g. realign,coreg,smooth): "
            read -r stages
            ;;
        [Xx])
            print_info "Cancelled."
            return 0
            ;;
        *)
            print_error "Invalid choice."
            return 1
            ;;
    esac

    stages="${stages// /}"
    if [[ -z "${stages}" ]]; then
        print_error "No stages given."
        return 1
    fi

    # Reject typos before anything reaches the cluster
    local valid=" topup vdm realign coreg smooth all fieldmap post_fieldmap "
    local part
    local IFS=','
    for part in ${stages}; do
        if [[ ! "${valid}" == *" ${part} "* ]]; then
            print_error "Unknown stage: '${part}'"
            printf "  Valid stages: topup, vdm, realign, coreg, smooth\n"
            printf "  Shorthands:   all, fieldmap, post_fieldmap\n"
            return 1
        fi
    done
    unset IFS

    if [[ "${PREPROC_MODE:-}" != "topup" ]] && [[ ",${stages}," == *",topup,"* ]]; then
        print_warning "PREPROC_MODE is '${PREPROC_MODE:-<not set>}', so the topup stage will do nothing."
    fi

    print_info "Stages selected: ${stages}"
    printf "\n"

    PREPROC_STAGES_OVERRIDE="${stages}"
    local rc=0
    run_step_2 || rc=$?
    PREPROC_STAGES_OVERRIDE=""
    return "${rc}"
}

# Step 3: Prepare data for GridCAT
run_step_3() {
    print_header "STEP 3: Prepare Data for GridCAT"

    if [[ ! -f "${SCRIPT_DIR}/run_prep.sh" ]]; then
        print_error "run_prep.sh not found in ${SCRIPT_DIR}"
        return 1
    fi

    print_info "This step will reorganize preprocessed data for GridCAT analysis."
    printf "  Script: run_prep.sh (which calls prepare_gridcat_directory.m)\n"
    printf "  Input:  Preprocessed SPM data\n"
    printf "  Output: GridCAT-ready directory structure\n\n"
    printf "NOTE: Make sure Step 2 preprocessing jobs have completed before running this step.\n"
    printf "You can check job status using the 'Check job status' menu option.\n\n"

    if ! ask_confirmation "Submit Step 3?"; then
        print_info "Step 3 skipped."
        return 0
    fi

    printf "\n"
    print_info "Submitting preparation job to cluster..."

    local -a dep_flag=()
    [[ -n "${NEXT_DEPENDENCY}" ]] && dep_flag=(--dependency="${NEXT_DEPENDENCY}")
    LAST_JOB_ID=""

    local sbatch_output
    if sbatch_output=$(sbatch \
        --export=ALL,PIPELINE_SCRIPT_DIR="${SCRIPT_DIR}" \
        --partition="${SLURM_PARTITION}" \
        --time="${PREP_TIME}" \
        --cpus-per-task="${PREP_CPUS}" \
        --mem="${PREP_MEM}" \
        "${dep_flag[@]}" \
        "${SCRIPT_DIR}/run_prep.sh" 2>&1); then
        print_success "Preparation job submitted!"

        local job_id
        job_id=$(echo "${sbatch_output}" | grep -oP 'Submitted batch job \K[0-9]+' || true)
        LAST_JOB_ID="${job_id}"

        if [[ -n "${job_id}" ]]; then
            printf "\n${GREEN}Job ID: ${job_id}${NC}\n"
            printf "\nYou can check the status with:\n"
            printf "  squeue -j ${job_id}\n\n"
            log_message "INFO" "Step 3 preparation job submitted with Job ID: ${job_id}"
        else
            printf "\n${GREEN}${sbatch_output}${NC}\n"
            log_message "INFO" "Step 3 preparation job submitted"
        fi

        printf "Please wait for this job to complete before running Step 4.\n"
        return 0
    else
        print_error "Failed to submit Step 3 job"
        printf "Error output:\n${sbatch_output}\n"
        log_message "ERROR" "Step 3 failed to submit job: ${sbatch_output}"
        return 1
    fi
}

# Step 4: Run GridCAT analysis
run_step_4() {
    print_header "STEP 4: Run GridCAT Analysis (SLURM)"

    if [[ ! -f "${SCRIPT_DIR}/run_gridcat.sh" ]]; then
        print_error "run_gridcat.sh not found in ${SCRIPT_DIR}"
        return 1
    fi

    print_info "This step will run the GridCAT GLM statistical analysis."
    printf "  Script: run_gridcat.sh (which calls run_gridcat_analysis.m)\n"
    printf "  Analysis type: Group-level GLM statistics\n"
    if [[ -n "${RUN_VARIANT:-}" ]]; then
        printf "  Run variant: ${RUN_VARIANT}  → output: ${OUTPUT_ROOT}/GLM_output_${RUN_VARIANT}\n"
    else
        printf "  Run variant: <none>  → output: ${OUTPUT_ROOT}/GLM_output\n"
    fi
    printf "\nNOTE: Make sure Step 3 has completed before running this step.\n"
    printf "NOTE: Each submission snapshots pipeline_config.cfg, so you can edit\n"
    printf "      RUN_VARIANT (and other settings) and submit another job while\n"
    printf "      this one is still running.\n\n"

    if ! ask_confirmation "Submit Step 4?"; then
        print_info "Step 4 skipped."
        return 0
    fi

    printf "\n"
    print_info "Submitting GridCAT analysis job to cluster..."

    # Delegate to run_gridcat.sh's submitter mode — it handles snapshotting
    # and builds the sbatch command from the config. NEXT_DEPENDENCY (if set by
    # run-all) is passed through GRIDCAT_DEPENDENCY so the GLM job waits for prep.
    LAST_JOB_ID=""
    local sbatch_output
    if sbatch_output=$(GRIDCAT_DEPENDENCY="${NEXT_DEPENDENCY}" bash "${SCRIPT_DIR}/run_gridcat.sh" 2>&1); then
        print_success "GridCAT analysis job submitted!"
        printf "%s\n" "${sbatch_output}"

        local job_id
        job_id=$(echo "${sbatch_output}" | grep -oP 'Submitted batch job \K[0-9]+' || true)
        LAST_JOB_ID="${job_id}"

        if [[ -n "${job_id}" ]]; then
            printf "\n${GREEN}Job ID: ${job_id}${NC}\n"
            printf "\nYou can check the status with:\n"
            printf "  squeue -j ${job_id}\n\n"
            log_message "INFO" "Step 4 GridCAT analysis job submitted with Job ID: ${job_id}"
        else
            printf "\n${GREEN}${sbatch_output}${NC}\n"
            log_message "INFO" "Step 4 GridCAT analysis job submitted"
        fi

        printf "GridCAT analysis is the final pipeline stage.\n"
        printf "Results will appear in: ${OUTPUT_ROOT}/GLM_output${RUN_VARIANT:+_${RUN_VARIANT}}/\n"
        return 0
    else
        print_error "Failed to submit Step 4 job"
        printf "Error output:\n${sbatch_output}\n"
        log_message "ERROR" "Step 4 failed to submit job: ${sbatch_output}"
        return 1
    fi
}

################################################################################
# DISPLAY FUNCTIONS
################################################################################

# Display the main menu
show_main_menu() {
    print_header "GridCAT Analysis Pipeline — Main Menu"

    printf "Individual Steps:\n"
    printf "  ${BLUE}[0]${NC}  Discover subjects & sessions\n"
    printf "  ${BLUE}[F]${NC}  Filter subject/session list\n"
    printf "  ${BLUE}[1]${NC}  Validate configuration\n"
    printf "  ${BLUE}[2]${NC}  Run SPM preprocessing (SLURM)\n"
    printf "  ${BLUE}[2S]${NC} Run SPM preprocessing — pick stages (SLURM)\n"
    printf "  ${BLUE}[3]${NC}  Prepare data for GridCAT (SLURM)\n"
    printf "  ${BLUE}[4]${NC}  Run GridCAT analysis (SLURM)\n\n"

    printf "Data Preparation:\n"
    printf "  ${BLUE}[M]${NC}  Scan for multi-run cases (T2w/BOLD)\n"
    printf "  ${BLUE}[R]${NC}  Move ROI masks into BIDS (dry-run)\n"
    printf "  ${BLUE}[RX]${NC} Move ROI masks into BIDS (execute)\n\n"

    printf "Batch Operations:\n"
    printf "  ${BLUE}[A]${NC}  Run ALL steps (0→1 locally, then 2→3→4 chained on SLURM)\n"
    printf "  ${BLUE}[S]${NC}  Show current settings\n"
    printf "  ${BLUE}[C]${NC}  Check job status\n"
    printf "  ${BLUE}[L]${NC}  Show pipeline log\n"
    printf "  ${BLUE}[Q]${NC}  Quit\n\n"

    printf "Enter your choice: "
}

# Show current settings from config
show_settings() {
    print_header "Current Pipeline Settings"

    if ! source_config; then
        print_error "Could not load configuration file"
        return 1
    fi

    printf "${BLUE}Data Paths:${NC}\n"
    printf "  BIDS Directory:        ${BIDS_ROOT:-<not set>}\n"
    printf "  Subject/Session List:  ${SCRIPT_DIR}/subses_list.txt\n"
    printf "  Output Directory:      ${OUTPUT_ROOT:-<not set>}\n"
    printf "  SPM Directory:         ${SPM_DIR:-<not set>}\n"
    printf "  GridCAT Directory:     ${GRIDCAT_DIR:-<not set>}\n"
    printf "  CircStat Directory:    ${CIRCSTAT_DIR:-<auto-detect>}\n\n"

    printf "${BLUE}Acquisition Parameters:${NC}\n"
    printf "  TR (sec):              ${TR:-<not set>}\n"
    printf "  Total Readout (ms):    ${TOTAL_READOUT_MS:-<not set>}\n"
    printf "  TE Short (ms):         ${TE_SHORT_MS:-<not set>}\n"
    printf "  TE Long (ms):          ${TE_LONG_MS:-<not set>}\n"
    printf "  Blip Direction:        ${BLIP_DIRECTION:-<not set>}\n"
    printf "  EPI-based Fieldmap:    ${EPI_BASED_FIELDMAP:-<not set>}\n\n"

    printf "${BLUE}Preprocessing Settings:${NC}\n"
    printf "  Preprocessing Mode:    ${PREPROC_MODE:-<not set>}\n"
    printf "  Preprocessing Stages:  ${PREPROC_STAGES:-all}\n"
    printf "  Coregister ROIs:       ${DO_COREG_ROIS:-<not set>}\n"
    printf "  Reslice Prefix:        ${RESLICE_PREFIX:-<not set>}\n\n"

    if [[ "${PREPROC_MODE:-}" == "topup" ]]; then
        printf "${BLUE}FSL Topup Settings:${NC}\n"
        printf "  FSL Module:            ${FSL_MODULE:-<not set>}\n"
        printf "  PE Dir (BOLD):         ${TOPUP_PE_DIR_BOLD:-<not set>}\n"
        printf "  PE Dir (Reverse):      ${TOPUP_PE_DIR_REVERSE:-<not set>}\n"
        printf "  Readout Time (sec):    ${TOPUP_READOUT_SEC:-<not set>}\n"
        printf "  Reverse-PE Search:     ${TOPUP_REVERSE_PE_DIR:-auto}\n"
        printf "  Reverse-PE Pattern:    ${TOPUP_REVERSE_PE_PATTERN:-<not set>}\n"
        printf "  Topup Config:          ${TOPUP_CONFIG:-<FSL default>}\n"
        printf "  Apply Method:          ${TOPUP_APPLY_METHOD:-<not set>}\n"
        printf "  Interpolation:         ${TOPUP_INTERP:-<not set>}\n"
        printf "  Existing topup prefix: ${TOPUP_EXISTING_PREFIX:-<use func/topup_results>}\n\n"
    fi

    printf "${BLUE}Study Design:${NC}\n"
    printf "  Tasks:                 ${TASKS:-<not set>}\n"
    printf "  Excluded Tasks:        ${EXCLUDE_TASKS:-<none>}\n"
    printf "  ROI Mode:              ${ROI_MODE:-<not set>}\n"
    printf "  Fold Symmetry:         ${X_FOLD_SYMMETRY:-<not set>}\n\n"

    printf "${BLUE}Subject/Session Filters:${NC}\n"
    printf "  Include Subjects:      ${INCLUDE_SUBJECTS:-<all>}\n"
    printf "  Include Sessions:      ${INCLUDE_SESSIONS:-<all>}\n"
    printf "  Exclude Subjects:      ${EXCLUDE_SUBJECTS:-<none>}\n"
    printf "  Exclude Sessions:      ${EXCLUDE_SESSIONS:-<none>}\n"
    printf "  Require ROI:           ${REQUIRE_ROI:-false}\n\n"

    printf "${BLUE}ROI & Multi-Run Settings:${NC}\n"
    printf "  ROI Source Dir:        ${ROI_SOURCE_DIR:-<not set>}\n"
    printf "  ROI Source Label:      ${ROI_SOURCE_LABEL:-<not set>}\n"
    printf "  ROI Pattern Left:      ${ROI_PATTERN_LEFT:-<not set>}\n"
    printf "  ROI Pattern Right:     ${ROI_PATTERN_RIGHT:-<not set>}\n"
    printf "  Run Selection File:    ${RUN_SELECTION_FILE:-<auto-select>}\n"
    if [[ -n "${RUN_SELECTION_FILE:-}" ]] && [[ -f "${RUN_SELECTION_FILE:-}" ]]; then
        local sel_count
        sel_count=$(tail -n +2 "$RUN_SELECTION_FILE" | grep -cve '^\s*$' || echo 0)
        printf "  Run Selection Entries: $sel_count\n"
    fi
    printf "\n"

    printf "${BLUE}GridCAT GLM Settings:${NC}\n"
    printf "  Masking Threshold:     ${MASKING_THRESHOLD:-<not set>}\n"
    printf "  HPF Cutoff (sec):      ${HPF_CUTOFF:-<not set>}\n"
    printf "  HRF Derivatives:       ${DERIVATIVES:-<not set>}\n"
    printf "  Func Prefix:           ${FUNC_PREFIX:-<not set>}\n"
    printf "  ROI Prefix:            ${ROI_PREFIX:-<not set>}\n\n"

    printf "${BLUE}HPC/SLURM Settings:${NC}\n"
    printf "  Partition:             ${SLURM_PARTITION:-<not set>}\n"
    printf "  MATLAB Module:         ${MATLAB_MODULE:-<not set>}\n"
    printf "  Run-all Chain Dep:     ${CHAIN_DEPENDENCY:-afterok}\n"
    printf "  Preproc Time/CPU/Mem:  ${PREPROC_TIME:-?} / ${PREPROC_CPUS:-?} CPUs / ${PREPROC_MEM_PER_CPU:-?} per CPU\n"
    printf "  Prep Time/CPU/Mem:     ${PREP_TIME:-?} / ${PREP_CPUS:-?} CPUs / ${PREP_MEM:-?}\n"
    printf "  GridCAT Time/CPU/Mem:  ${GRIDCAT_TIME:-?} / ${GRIDCAT_CPUS:-?} CPUs / ${GRIDCAT_MEM:-?}\n"
    printf "  Max Workers:           ${MAX_WORKERS:-<not set>}\n"
    printf "  Local scratch:         ${USE_LOCAL_SCRATCH:-false}${GRIDCAT_LOCAL_TMP:+ (${GRIDCAT_LOCAL_TMP})}\n\n"

    printf "${BLUE}Pipeline Behavior:${NC}\n"
    printf "  Copy Mode:             ${COPY_MODE:-<not set>}\n"
    printf "  Fail on Missing:       ${FAIL_ON_MISSING:-<not set>}\n"
    printf "  Dry Run:               ${DRY_RUN:-<not set>}\n\n"

    printf "${BLUE}Pipeline Directories:${NC}\n"
    printf "  Script Directory:      ${SCRIPT_DIR}\n"
    printf "  Config File:           ${CONFIG_FILE}\n"
    printf "  Log Directory:         ${LOG_DIR}\n"
    printf "  Current Log:           ${LOG_FILE}\n\n"
}

# Check job status with squeue
check_job_status() {
    print_header "HPC Job Status"

    printf "Checking running jobs for user: ${USER}\n\n"

    if squeue -u "${USER}"; then
        printf "\n"
        return 0
    else
        print_warning "Could not query job status. squeue may not be available."
        printf "This is normal if you're not on an HPC cluster.\n"
        return 1
    fi
}

# Show the pipeline log
show_pipeline_log() {
    print_header "Pipeline Activity Log"

    if [[ ! -f "${LOG_FILE}" ]]; then
        print_warning "No log file exists yet."
        return 0
    fi

    printf "Log file: ${LOG_FILE}\n\n"
    tail -50 "${LOG_FILE}"
    printf "\n"
}

# Submit the cluster stages (2 → 3 → 4) as a SLURM dependency chain.
# Local stages (0, filter, 1) run first; validation is a hard gate. Each
# submit step records its job id in LAST_JOB_ID, which becomes the next
# stage's --dependency. Returns non-zero (and stops submitting) on any failure.
_run_all_chain() {
    local chain_dep="${CHAIN_DEPENDENCY:-afterok}"
    if [[ "${chain_dep}" != "afterok" && "${chain_dep}" != "afterany" ]]; then
        print_warning "CHAIN_DEPENDENCY='${chain_dep}' is unusual; expected afterok or afterany."
    fi

    # ---- Local stages: run now ----
    if ! run_step_0; then
        print_error "Run-all aborted: Step 0 (subject discovery) failed."
        return 1
    fi
    printf "\n"

    run_step_filter
    printf "\n"

    if ! run_step_1; then
        print_error "Run-all aborted: Step 1 (validation) failed."
        printf "No cluster jobs were submitted. Fix the issues above and retry.\n"
        return 1
    fi
    printf "\n"

    # ---- Cluster stages: queue as a dependency chain ----
    NEXT_DEPENDENCY=""
    if ! run_step_2 || [[ -z "${LAST_JOB_ID}" ]]; then
        print_error "Run-all aborted: Step 2 did not submit (no job ID captured)."
        return 1
    fi
    local preproc_id="${LAST_JOB_ID}"
    printf "\n"

    NEXT_DEPENDENCY="${chain_dep}:${preproc_id}"
    if ! run_step_3 || [[ -z "${LAST_JOB_ID}" ]]; then
        print_error "Run-all aborted: Step 3 did not submit."
        printf "Step 2 (job ${preproc_id}) is still queued. Cancel with: scancel ${preproc_id}\n"
        return 1
    fi
    local prep_id="${LAST_JOB_ID}"
    printf "\n"

    NEXT_DEPENDENCY="${chain_dep}:${prep_id}"
    if ! run_step_4 || [[ -z "${LAST_JOB_ID}" ]]; then
        print_error "Run-all aborted: Step 4 did not submit."
        printf "Jobs ${preproc_id} and ${prep_id} are still queued. Cancel with: scancel ${preproc_id} ${prep_id}\n"
        return 1
    fi
    local gridcat_id="${LAST_JOB_ID}"
    printf "\n"

    # ---- Summary ----
    print_success "Full pipeline queued as a SLURM dependency chain."
    printf "\n"
    printf "  ${GREEN}Stage 2 — preprocessing : job %s${NC}\n" "${preproc_id}"
    printf "  ${GREEN}Stage 3 — prep          : job %s   (waits %s:%s)${NC}\n" "${prep_id}" "${chain_dep}" "${preproc_id}"
    printf "  ${GREEN}Stage 4 — GridCAT GLM   : job %s   (waits %s:%s)${NC}\n" "${gridcat_id}" "${chain_dep}" "${prep_id}"
    printf "\nSLURM starts each stage automatically once the previous one finishes —\n"
    printf "you do not need to come back and launch the next stage by hand.\n\n"
    printf "Monitor:    squeue -u ${USER}    (or menu option C)\n"
    printf "Cancel all: scancel ${preproc_id} ${prep_id} ${gridcat_id}\n"
    log_message "INFO" "Run-all chain submitted: preproc=${preproc_id} prep=${prep_id} gridcat=${gridcat_id} (dep=${chain_dep})"
    return 0
}

# Run all steps: confirm once, then run the chain non-interactively.
run_all_steps() {
    print_header "Run ALL Steps — chained SLURM submission"

    local chain_dep="${CHAIN_DEPENDENCY:-afterok}"
    printf "Local stages run now; cluster stages are queued as one dependency\n"
    printf "chain, so the whole pipeline runs unattended:\n\n"
    printf "  0  Discover subjects/sessions     ${BLUE}(now)${NC}\n"
    printf "  F  Apply filters                  ${BLUE}(now)${NC}\n"
    printf "  1  Validate configuration         ${BLUE}(now, hard gate)${NC}\n"
    printf "  2  SPM preprocessing   ${BLUE}─┐${NC}\n"
    printf "  3  Prepare for GridCAT ${BLUE} ├─ chained via --dependency=%s${NC}\n" "${chain_dep}"
    printf "  4  GridCAT GLM         ${BLUE}─┘${NC}\n\n"
    printf "If validation fails, no cluster jobs are submitted.\n\n"

    if ! ask_confirmation "Submit the full chain?"; then
        print_info "Run-all cancelled."
        return 0
    fi
    printf "\n"

    # Auto-approve the per-step prompts (we just confirmed the whole run),
    # then always restore interactive state regardless of outcome.
    ASSUME_YES=true
    local rc=0
    _run_all_chain || rc=$?
    ASSUME_YES=false
    NEXT_DEPENDENCY=""
    LAST_JOB_ID=""
    return "${rc}"
}

################################################################################
# MAIN MENU LOOP
################################################################################

main() {
    # Check if config file exists before starting
    if ! check_config; then
        exit 1
    fi

    # Try to source config for initial setup
    if ! source_config; then
        print_warning "Could not load configuration file. Some features may not work."
    fi

    log_message "INFO" "Pipeline runner started by user: ${USER}"

    local choice

    while true; do
        show_main_menu
        read -r choice

        case "${choice}" in
            0)
                run_step_0
                printf "\nPress Enter to continue..."
                read -r
                ;;
            [Ff])
                run_step_filter
                printf "\nPress Enter to continue..."
                read -r
                ;;
            1)
                run_step_1
                printf "\nPress Enter to continue..."
                read -r
                ;;
            [2][Ss])
                run_step_2_stages
                printf "\nPress Enter to continue..."
                read -r
                ;;
            2)
                run_step_2
                printf "\nPress Enter to continue..."
                read -r
                ;;
            3)
                run_step_3
                printf "\nPress Enter to continue..."
                read -r
                ;;
            4)
                run_step_4
                printf "\nPress Enter to continue..."
                read -r
                ;;
            [Aa])
                run_all_steps
                printf "\nPress Enter to continue..."
                read -r
                ;;
            [Ss])
                show_settings
                printf "Press Enter to continue..."
                read -r
                ;;
            [Cc])
                check_job_status
                printf "Press Enter to continue..."
                read -r
                ;;
            [Ll])
                show_pipeline_log
                printf "Press Enter to continue..."
                read -r
                ;;
            [Mm])
                print_header "Scanning for Multi-Run Cases"
                log_message "INFO" "Running multi-run scanner"
                bash "${SCRIPT_DIR}/scan_multirun.sh" "${CONFIG_FILE}" 2>&1 | tee -a "$LOG_FILE"
                printf "\nPress Enter to continue..."
                read -r
                ;;
            [Rr][Xx])
                print_header "Move ROI Masks (EXECUTE)"
                printf "${YELLOW}This will copy and rename ROI files into BIDS anat/ folders.${NC}\n"
                printf "Are you sure? (y/N): "
                read -r confirm
                if [[ "${confirm}" =~ ^[Yy]$ ]]; then
                    log_message "INFO" "Running ROI mover (execute)"
                    bash "${SCRIPT_DIR}/move_rois.sh" --execute 2>&1 | tee -a "$LOG_FILE"
                else
                    print_info "Cancelled."
                fi
                printf "\nPress Enter to continue..."
                read -r
                ;;
            [Rr])
                print_header "Move ROI Masks (dry-run)"
                log_message "INFO" "Running ROI mover (dry-run)"
                bash "${SCRIPT_DIR}/move_rois.sh" 2>&1 | tee -a "$LOG_FILE"
                printf "\nPress Enter to continue..."
                read -r
                ;;
            [Qq])
                print_info "Exiting pipeline runner."
                log_message "INFO" "Pipeline runner exited by user"
                exit 0
                ;;
            *)
                print_error "Invalid choice. Please enter 0-4, 2S, F, A, S, C, L, M, R, RX, or Q."
                printf "\nPress Enter to continue..."
                read -r
                ;;
        esac
    done
}

# Run main function
main "$@"
