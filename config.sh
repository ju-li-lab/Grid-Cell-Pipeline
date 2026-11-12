# =============================================================================
#  Grid Cell Pipeline — shell settings
# =============================================================================
#  Paths, cluster resources and which subjects/runs to process. Sourced by the
#  shell scripts, so plain bash: KEY=value, no spaces around '='.
#  SPM and GridCAT parameters live in spm_options.m, not here.
# =============================================================================


# -----------------------------------------------------------------------------
#  Paths
# -----------------------------------------------------------------------------

# BIDS data. Read only — the pipeline never writes here.
BIDS_ROOT=/sc-projects/sc-proj-cc02-brainspace/myspace/julius/SOLIS_GC/bids

# Where output goes. Created if missing.
OUTPUT_ROOT=/sc-projects/sc-proj-cc02-brainspace/myspace/julius/SOLIS_GC/workDir

# This script directory on the cluster. run_selection.tsv lives here.
PIPELINE_DIR=/sc-projects/sc-proj-cc02-brainspace/myspace/julius/SOLIS_GC/Grid-Cell-Pipeline

SPM_DIR=/sc-projects/sc-proj-cc02-brainspace/software/matlab/spm25
GRIDCAT_DIR=/sc-projects/sc-proj-cc02-brainspace/myspace/julius/sirius_wip/GridCAT

# Ships with GridCAT. Empty = GRIDCAT_DIR/CircStat2012a
CIRCSTAT_DIR=


# -----------------------------------------------------------------------------
#  Derivatives — everything the pipeline writes
# -----------------------------------------------------------------------------
#  SPM writes estimated transforms into the headers of its INPUT images, so
#  every input is copied here first and the copy is what gets touched.

# Empty = <OUTPUT_ROOT>/derivatives
DERIV_ROOT=

# Three datasets inside it. spm-preproc is regenerated on every full run
# (see DERIV_RESET); rois and fieldmaps are kept.
DERIV_PREPROC=spm-preproc
DERIV_ROIS=rois
DERIV_FIELDMAPS=fieldmaps

# When to clear the preprocessing derivatives:
#   auto    clear when the run starts from its first stage; a resume keeps them
#   always  clear on every run (resuming stops working)
#   never   never clear; new output overwrites old file by file
DERIV_RESET=auto

# Names this run's GridCAT output folder: empty -> GLM_output, otherwise
# GLM_output_<variant>. Letters, digits, underscore, dash and dot only.
# The settings behind each variant are recorded in its run_info.tsv.
RUN_VARIANT=6f_t5_am_keep1_avg_1


# -----------------------------------------------------------------------------
#  What to process
# -----------------------------------------------------------------------------

# Task labels, comma-separated, matching task-<LABEL> in the BIDS filenames.
TASKS=1,2,3

# How to correct distortion:
#   realign_unwarp   GRE fieldmap: needs phasediff + magnitude1 in fmap/
#   topup            FSL topup from a reverse-PE EPI
#   precalc_fieldmap ready-made B0 fieldmap (tools/make_fieldmaps.sh)
#   realign_only     motion correction only
PREPROC_MODE=realign_unwarp

# Which preprocessing stages to run, comma-separated: topup, vdm, slicetime,
# realign, coreg, smooth. Each reads what it needs from disk, so you can stop
# after any stage and resume later. Stages always run in pipeline order.
# Shorthands: all, fieldmap (= topup,vdm), post_fieldmap (= everything after).
# Override for one submission:  bash submit_preproc.sh --stages smooth
PREPROC_STAGES=all

# Coregister T2w -> EPI and reslice ROI masks into EPI space? (1=yes, 0=skip)
DO_COREG_ROIS=1

# Coregister the T2w to the session mean THROUGH a whole-brain EPI with this
# task label ('rest' finds func/*_task-rest*_bold.nii). Empty = direct route.
# Use it when T2w and task EPI are both partial slabs: their small overlap can
# make NMI diverge (sub-05 ses-01 ended 141 mm off with empty masks). The
# masks are still resliced only once, with the composed transform.
COREG_VIA_TASK=rest

# Fail the job when a resliced ROI mask has fewer voxels than this. A diverged
# coregistration writes a mask that misses the brain instead of erroring.
# 0 disables the check.
COREG_MIN_ROI_VOXELS=3


# -----------------------------------------------------------------------------
#  Subject / session filtering  (00_list_subjects.sh; empty = everything)
# -----------------------------------------------------------------------------

INCLUDE_SUBJECTS=
INCLUDE_SESSIONS=
EXCLUDE_SUBJECTS=sub-07s17
EXCLUDE_SESSIONS=

# Skip subject-sessions that have no ROI masks in anat/.
REQUIRE_ROI=true

# Tasks to exclude from the GridCAT analysis, comma-separated.
EXCLUDE_TASKS=rest,reverse,mean


# -----------------------------------------------------------------------------
#  ROI masks
# -----------------------------------------------------------------------------

# Filenames the pipeline looks for in anat/. Extra BIDS entities tolerated.
ROI_PATTERN_LEFT=_hemi-left_label-lpmErC_mask.nii
ROI_PATTERN_RIGHT=_hemi-right_label-rpmErC_mask.nii

# The anatomical label in those filenames. Everything downstream filters on
# this, so it must match the patterns.
ROI_LABEL=ErC

# Where tools/move_rois.sh imports manual segmentations from. Files are
# expected as: sub-XX_ses-YY_[run-N_]<ROI_LABEL>_<hemi>.nii[.gz]
ROI_SOURCE_DIR=/sc-projects/sc-proj-cc02-brainspace/myspace/julius/BRACE_GC_testing/v2_pipe_t5_data/ROI_lucie/ErC_ROIs_BIDS

# ROI analysis mode for GridCAT GLM2: both | bilat_only | lr_only
ROI_MODE=both

# Use tSNR-filtered ROI masks for the GridCAT GLM2? (1=yes, 0=no)
# 1 makes 03_gridcat.sh read masks from GLM_runauto/ROI_tsnr, written by
# tools/roi_tsnr/run_roi_tsnr_filter.sh (run that first). 0 uses GLM_runauto/ROI.
ROI_TSNR_FILTER=0


# -----------------------------------------------------------------------------
#  Multiple runs of the same scan
# -----------------------------------------------------------------------------
#  Run numbers say nothing across modalities, so the pipeline decides:
#    1. your choice in <PIPELINE_DIR>/run_selection.tsv (bash tools/scan_multirun.sh)
#    2. the run acquired closest in time (AcquisitionTime in the JSONs)
#    3. nothing else is safe — see STRICT_RUN_MATCHING

# How far apart two scans may be acquired and still count as one visit, minutes.
RUN_MATCH_GAP_MIN=20

# When nothing can decide between several runs:
#   true   stop and say so (default — a wrong pairing is invisible afterwards)
#   false  take the highest run number and warn loudly
STRICT_RUN_MATCHING=true


# -----------------------------------------------------------------------------
#  Finding fieldmap inputs
# -----------------------------------------------------------------------------
#  Filename patterns only. The numbers that go into SPM are in spm_options.m.

# FSL module name, for: module load <this>
FSL_MODULE=fsl/6.0.7.18

# --- topup ---
# Where the reverse-PE EPI is: fmap | func | auto (fmap first, then func)
TOPUP_REVERSE_PE_DIR=func

# Filename substring for the reverse-PE EPI.
TOPUP_REVERSE_PE_PATTERN=task-reverse_bold

# FSL topup config file. Empty = FSL default b02b0.cnf
TOPUP_CONFIG=

# How to apply the correction:
#   vdm         convert the topup field to an SPM VDM, then Realign & Unwarp
#   applytopup  correct the EPIs with FSL applytopup, then SPM Realign only
TOPUP_APPLY_METHOD=vdm

# Pre-computed topup output, used when the "topup" stage is NOT run. Prefix
# without the _fieldcoef.nii.gz suffix; {SUB}/{SES} substituted. Empty = use
# func/topup_results.
TOPUP_EXISTING_PREFIX=

# --- precalc_fieldmap ---
# Patterns for fmap/, WITHOUT the extension.
FIELDMAP_PATTERN=_fieldmap
FIELDMAP_MAGNITUDE_PATTERN=_magnitude

# Units on disk: rad/s (fsl_prepare_fieldmap output) or Hz. rad/s is divided
# by 2*pi before use — getting this wrong scales the correction by 2*pi.
FIELDMAP_UNITS=rad/s


# -----------------------------------------------------------------------------
#  Building fieldmaps from a structural (tools/make_fieldmaps.sh)
# -----------------------------------------------------------------------------
#  For sessions with a phasediff but no magnitude: builds the magnitude from
#  the T1w. Run with --dry-run first, then set PREPROC_MODE=precalc_fieldmap.

# Passed straight to fsl_prepare_fieldmap: SIEMENS or PHILIPS
FIELDMAP_SCANNER=SIEMENS

# Echo time difference in ms. Empty = read per session from the phasediff JSON.
FIELDMAP_DELTA_TE=

# bet fractional intensity threshold (lower keeps more brain).
FIELDMAP_BET_F=0.5

# Where the magnitude comes from:
#   structural  always build from the T1w (every session processed identically)
#   auto        use the real *_magnitude1 when there is one, T1w otherwise
FIELDMAP_MAGNITUDE_SOURCE=structural


# -----------------------------------------------------------------------------
#  Output file prefixes — must match what SPM writes
# -----------------------------------------------------------------------------
#  The prefix of the final BOLD series is derived from these plus the switches
#  below (func_prefix() in bids_common.sh), so nothing to keep in sync by hand.

RESLICE_PREFIX=u
SMOOTH_PREFIX=s
ROI_PREFIX=r
SLICETIME_PREFIX=a

# Slice timing correction? (1=yes, 0=no). Reads SliceTiming from each run's
# JSON sidecar and fails loudly if absent. Changing this changes which files
# later stages read — re-run from "slicetime" onwards.
DO_SLICE_TIMING=1

# Slice timing before or after realignment: before -> ua*, after -> au*.
# Both are correct SPM usage; see README. Flipping it redoes both operations
# (--stages post_fieldmap), and the two orders belong in separate trees.
SLICETIME_ORDER=before

# Smoothing kernel FWHM in mm, isotropic. 0 skips smoothing entirely.
# Lives here because it decides which files later stages read (su* vs u*).
SMOOTH_FWHM=0


# -----------------------------------------------------------------------------
#  SLURM
# -----------------------------------------------------------------------------

SLURM_PARTITION=compute

# Job notification email. Empty disables it.
NOTIFY_EMAIL=julius.schaut@charite.de

MATLAB_MODULE=matlab

# How run_all.sh chains the jobs: afterok (next stage only if all succeeded)
# or afterany.
CHAIN_DEPENDENCY=afterok

# Preprocessing: one array task per subject-session
PREPROC_TIME=10:00:00
PREPROC_CPUS=2
PREPROC_MEM_PER_CPU=8G

# GridCAT preparation
PREP_TIME=01:00:00
PREP_CPUS=10
PREP_MEM=40G

# GridCAT analysis
GRIDCAT_TIME=01:00:00
GRIDCAT_CPUS=30
GRIDCAT_MEM=80G

# Parallel MATLAB workers inside the GridCAT job, capped at the number of
# subject-sessions. Rule of thumb: GRIDCAT_CPUS minus 1.
MAX_WORKERS=29

# Stage data to node-local scratch to bypass NFS?
USE_LOCAL_SCRATCH=true

# Local tmp to request from SLURM (--gres=tmp:SIZE). Must exceed GLM_runauto
# plus expected output. Empty skips the request.
GRIDCAT_LOCAL_TMP=80G

# Where to stage on the node. Empty = SLURM $TMPDIR or /tmp.
LOCAL_SCRATCH_BASE=


# -----------------------------------------------------------------------------
#  Behaviour
# -----------------------------------------------------------------------------

# How files enter the GridCAT directory: copy (safer) | symlink (saves space)
COPY_MODE=copy

# On a missing file: true = stop so you can fix it, false = skip and continue
FAIL_ON_MISSING=false
