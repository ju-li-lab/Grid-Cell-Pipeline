+++ 09.07.2026 - UPDATE: Moved Pipeline to GitHub repository. To run the pipeline, clone using http link and run as usual. +++


Version ==13.03.2026== edited 27.04.2026

This is still a work in progress. Script paths will be provided.
For now you can make a copy of the folder into your own directory

```
/sc-projects/sc-proj-cc02-brainspace/sirius_wip/v2_fR/scripts_v2
```

! Usage of "session", "task", and "run" can differ between documentations. For example the "sessions" in GridCat don't match what is meant by "session" in bids but instead refers to "tasks"


---

## Overview

This document describes the complete GridCAT analysis pipeline, including SPM-based preprocessing, motion parameter extraction, ROI registration, and GLM-based connectivity analysis.

The pipeline is built around a **single configuration file** (`pipeline_config.cfg`) that controls all processing steps. A **master runner script** (`run_pipeline.sh`) provides an interactive menu to execute any or all pipeline steps in sequence.

All scripts auto-detect their directory and read configuration automatically.

---

## Quick Start

1. **Prepare your BIDS dataset** with anat (T2w + ROI masks), fmap, and func directories
2. **Copy the pipeline** to your project directory (do not modify the original scripts)
3. **Edit `pipeline_config.cfg`** with your project-specific settings (dataset path, task names, ROI labels, HPC options, etc.)
4. **Run the master pipeline**:
   ```bash
   ./run_pipeline.sh
   ```
5. **Select steps** from the interactive menu (0-5)
6. **Monitor progress** via console output and HPC logs
7. **Collect results** using the R output script (currently not working)

---

## Files in This Pipeline

| File                          | Type   | Purpose                                                                        |
| ----------------------------- | ------ | ------------------------------------------------------------------------------ |
| `run_pipeline.sh`             | Shell  | Master runner; interactive menu to execute pipeline steps                      |
| `pipeline_config.cfg`         | Config | Central configuration file; user edits this only                               |
| `step0_make_subses_list.sh`   | Shell  | Discovers all subject/session pairs from BIDS directory                        |
| `filter_subses_list.sh`       | Shell  | Filters subject/session list by include/exclude rules and ROI availability     |
| `validate_config.sh`          | Shell  | Validates all config settings before running analysis                          |
| `spm_preproc_array.sbatch`    | SLURM  | HPC array job script for SPM preprocessing (one per subject-session)           |
| `run_spm_preproc.m`           | MATLAB | Performs VDM calc, realign/unwarp (or realign-only), coregistration, reslicing |
| `read_pipeline_config.m`      | MATLAB | Config parser; reads and validates settings in MATLAB                          |
| `run_prep.sh`                 | Shell  | Wrapper for data preparation for GridCAT                                       |
| `prepare_gridcat_directory.m` | MATLAB | Splits 4D → 3D, copies motion regressors, creates bilateral ROIs               |
| `run_gridcat.sh`              | Shell  | Wrapper for GridCAT GLM analysis                                               |
| `run_gridcat_analysis.m`      | MATLAB | Runs GLM1 (per-session) and GLM2 (per-ROI) connectivity analysis               |
| `collect_gridcat_output.R`    | R      | Collects GridCAT results into CSV + HTML report                                |
| `move_rois.sh`                | Shell  | Copies/renames manual ROI segmentations into BIDS anat/ directories            |
| `scan_multirun.sh`            | Shell  | Scans for multi-run T2w/BOLD/reverse-PE files; generates `run_selection.tsv`   |
| `ashs_extract_config.sh`      | Shell  | (Optional) Helper to extract ASHS ROI configuration                            |
| `ashs_extract_labels.sh`      | Shell  | (Optional) Extracts ROI labels from ASHS output                                |

---

## Assumptions & Input Requirements

### Minimal BIDS Structure

```
dataset/
  sub-XX/
    ses-YY/
      anat/
        sub-XX_ses-YY_T2w.nii.gz          (anatomical reference)
        sub-XX_ses-YY_T2w_hemi-left_label-ErC_mas.nii.gz  (left ROI label-valued)
        sub-XX_ses-YY_T2w_hemi-right_label-ErC_mas.nii.gz (right, label-valued)
      fmap/
        sub-XX_ses-YY_magnitude1.nii.gz    (GRE FieldMap; for realign_unwarp)
        sub-XX_ses-YY_phasediff.nii.gz     (GRE FieldMap; for realign_unwarp)
        sub-XX_ses-YY_phasediff.json
      func/
        sub-XX_ses-YY_task-1_bold.nii.gz
        sub-XX_ses-YY_task-1_bold.json
        sub-XX_ses-YY_task-2_bold.nii.gz
        sub-XX_ses-YY_task-2_bold.json
        sub-XX_ses-YY_task-reverse_bold.nii.gz   (reverse-PE EPI; for topup mode)
        sub-XX_ses-YY_task-reverse_bold.json
```

**Note on fieldmap requirements by preprocessing mode:**
- `realign_unwarp`: Needs `magnitude1` + `phasediff` in fmap/
- `topup`: Needs a reverse phase-encode EPI (can live in fmap/ or func/, configurable via `TOPUP_REVERSE_PE_DIR` and `TOPUP_REVERSE_PE_PATTERN`)
- `realign_only`: No fieldmap needed

### Critical Assumptions

- **Each run has exactly one motion regressor file** (rp*.txt), or search pattern must be unique.
- **All tasks for the same GLM must have comparable acquisition parameters** (phase-encode direction, readout time, voxel size, TR).
- **ROI masks are label masks** (integer-valued), not probability maps.
- **You know the space** of every file at each processing stage (native T2w, native EPI, unwarped EPI, resliced EPI).
- **Runs with different acquisition parameters** require explicit handling (branching, separate preprocessing, or exclusion).
- **Multiple runs of the same task** are supported. Use `scan_multirun.sh` to generate a `run_selection.tsv` that lets you choose which run to use per subject/session. Without a selection file, the pipeline auto-selects the last (highest-numbered) run.

---

## Global Caveats

### Never Preprocess the Same File Twice

- Do not rerun realign/unwarp/coregister on a product of a previous preprocessing step.
- Always start from a fresh, native copy of the functional data.
- If uncertain, delete the derived directory and rebuild from raw/copies via bidscoiner or manual reset.
- Cumulative alterations can produce subtle artifacts. When in doubt, inspect headers and provenance.

### Make Study Design Explicit Upfront

Before running the pipeline, clarify for your project:

- **Which tasks are included in GridCAT?** Which are excluded (rest, other conditions)?
- **Unwarp, topup, or realign-only?** Set `PREPROC_MODE` in config. If topup, also set `TOPUP_APPLY_METHOD`.
- **Are all runs identical?** Different phase-encode directions, readout times, TR, multiband factors, voxel sizes?
- **Non-uniform runs** (different number of volumes, resting-state vs task, etc.)?

If runs are not homogeneous, scripts need manual review and possible branching logic.

---

## Extracting ROI masks from ASHS

If you already have ashs segmentations for your data, they will likely contain more than one label within the output file. The script directory contains two scripts that help you extract the labels you want to include in your GLM. Since the Labels have to be coregistered to the task space, you should do this first.

### Using the provided extraction scripts

For easier use, the script works with a config file that lets you specify all important paths and settings. There are some things you should know before running:

1. Specify where your ashs segmentations live -> The script assumes you are giving it the direct path to the atlas, so don't choose the general ashs folder if you ran it with more than one atlas.
2. You will have to start a slurm using something like ``` srun --pty bash ``` before running the actual script with ```bash ashs_extract_labels.sh ```
3. After confirming that the masks ended up in the right place, don't forget to change the ROI patterns in the general config file accordingly.

---

## Configuration: pipeline_config.cfg

The `pipeline_config.cfg` file controls all pipeline behavior. Edit this file once per project; all scripts read it automatically.

### Structure

The config file contains 7 main sections:

1. **PATHS** — BIDS_ROOT, OUTPUT_ROOT, SCRIPT_DIR, SPM_DIR, GRIDCAT_DIR
2. **ACQUISITION PARAMETERS** — Readout time, echo times, blip direction, TR
3. **STUDY DESIGN** — Tasks, preprocessing mode, ROI patterns, exclusions
4. **SPM PREPROCESSING SETTINGS** — Realign, unwarp, reslice, VDM, FSL topup, coregistration parameters
5. **GRIDCAT GLM SETTINGS** — Fold symmetry, masking, microtime, HPF, derivatives, event usage, contrast mode
6. **SLURM / HPC SETTINGS** — Partition, MATLAB module, time/CPU/memory per job step
7. **PIPELINE BEHAVIOR** — Copy mode, fail-on-missing, dry-run

### Key Settings to Review

- **BIDS_ROOT** — Absolute path to BIDS dataset root
- **OUTPUT_ROOT** — Where to write derived data
- **PREPROC_MODE** — Choose `realign_unwarp`, `topup`, or `realign_only`
- **TASKS** — Comma-separated task numbers to process (e.g., `1,2,3` matches task-1, task-2, task-3)
- **ROI_MODE** — Which ROI masks to use (`both`, `bilat_only`, or `lr_only`)
- **FAIL_ON_MISSING** — If `true`, pipeline stops on any missing file; if `false`, lenient processing
- **DRY_RUN** — If `true`, scripts preview actions without running (useful for validation)
- **SLURM_PARTITION** — HPC partition name (check with `sinfo -s`)

All preprocessing, GLM, and HPC parameters are configurable in `pipeline_config.cfg`. Review the advanced SPM settings (section 4 in the config) if you need to tune realignment, unwarping, or coregistration behavior for your specific data.

### Subject/Session Filtering

The pipeline supports flexible filtering of which subject-session pairs to process. This is configured in `pipeline_config.cfg` and applied by `filter_subses_list.sh` (run automatically after Step 0, or manually via the `[F]` menu option).

**Filter settings:**

| Setting | Type | Description |
|---------|------|-------------|
| `INCLUDE_SUBJECTS` | Whitelist | Only process these subjects (comma-separated, e.g., `sub-01s13,sub-02s13`) |
| `INCLUDE_SESSIONS` | Whitelist | Only process these sessions (comma-separated, e.g., `ses-01,ses-04`) |
| `EXCLUDE_SUBJECTS` | Blacklist | Remove these subjects (comma-separated) |
| `EXCLUDE_SESSIONS` | Blacklist | Remove these sessions (comma-separated) |
| `REQUIRE_ROI` | Boolean | If `true`, only keep subject-sessions that have ROI masks in `anat/` |

**How filters are applied:**
1. INCLUDE filters are applied first (whitelist)
2. EXCLUDE filters are applied second (blacklist removes from the whitelist result)
3. REQUIRE_ROI is applied last (removes pairs without ROI masks)

**Examples:**
```bash
# Only process sessions 1 and 4 for all subjects
INCLUDE_SESSIONS=ses-01,ses-04

# Process all subjects except sub-07s13
EXCLUDE_SUBJECTS=sub-07s13

# Only process subjects that have ROI masks
REQUIRE_ROI=true

# Combine: only specific subjects, only sessions with ROIs
INCLUDE_SUBJECTS=sub-01s13,sub-02s13,sub-03s13
REQUIRE_ROI=true
```

Leave all filter settings empty to process all subject-session pairs (default behavior).

The filter creates a backup of the original list as `subses_list_full.txt` before modifying `subses_list.txt`. The filtered list is then used by preprocessing (Step 2), GridCAT preparation (Step 3), and GridCAT analysis (Step 4).

---

## Step-by-Step Workflow

### Step 0: Discover Subjects and Sessions

```bash
./step0_make_subses_list.sh
```

**What it does:**
- Scans the BIDS dataset root for all `sub-XX/ses-YY/` directories
- Writes `subses_list.txt` in the pipeline root (one subject-session per line)
- Used by HPC to create the array job range

**Output:**
- `subses_list.txt` (plain text, format: `sub-XX_ses-YY`)

**QC check:**
- Verify the list contains all expected subjects/sessions
- Check for typos in naming or duplicate entries

---

### Option M: Scan for Multi-run cases (recommended)

If you have not pre-filtered and removed cases with more than one run of a scan, use this tool.

**What it does:**
- Creates a filter csv file 
- CSV file is used downstream for run selection

**How to use**
- Run the step
- Edit the csv file using ```nano file_name.csv```

**!CAVE!** Do NOT rely entirely upon this script. This is a messy situation and can be easy to disregard because of how much work it takes, but you absolutely need to manually double check this to make sure that:
- You pass uncorrupted and high quality images into the pipeline
- You consistently use either a corrected or uncorrected version of a scan across the study
- Your anatomical images match your functional scans I.E. the subject cannot have left the scanner between the aquisition of T2w and the functional scans!

---

### Step 0b: Filter Subject/Session List (Optional)

```bash
./filter_subses_list.sh
```

Or use the `[F]` option in the interactive menu. This step runs automatically when using "Run ALL steps".

**What it does:**
- Reads `subses_list.txt` (from Step 0)
- Backs up the original list as `subses_list_full.txt`
- Applies INCLUDE filters (whitelist), then EXCLUDE filters (blacklist), then REQUIRE_ROI
- Writes the filtered result back to `subses_list.txt`
- Prints a summary of how many pairs were kept vs removed

**When to use:**
- You only need to process specific sessions (e.g., `ses-01` and `ses-04`)
- You want to exclude certain subjects (e.g., bad data quality)
- You want to skip subject-sessions without ROI masks (saves preprocessing time)
- You have a large dataset (300+ combinations) and only need a subset

**Output:**
- `subses_list.txt` (filtered, used by all downstream steps)
- `subses_list_full.txt` (backup of original unfiltered list)

**Note:** If no filter settings are configured in `pipeline_config.cfg`, this step does nothing and all pairs are processed.

---

### Step 1: Validate Configuration

```bash
./validate_config.sh
```

**What it does:**
- Reads `pipeline_config.cfg`
- Checks that all required settings are defined
- Verifies that referenced directories/files exist
- Tests MATLAB and R availability (if needed)
- Checks SLURM/HPC connectivity (if HPC mode enabled)

**Output:**
- Console messages indicating pass/fail for each check
- If failures occur, script exits with error message

**When to run:**
- Before the first pipeline run
- After any major config change
- To diagnose "missing file" or "bad path" errors

---

### Step 2: SPM Preprocessing (HPC)

```bash
sbatch --array=1-N spm_preproc_array.sbatch
```

(Where `N` is the number of lines in `subses_list.txt`)

Or use the interactive menu in `run_pipeline.sh` to submit.

**What it does:**

Each job processes one subject-session pair. The exact steps depend on `PREPROC_MODE`:

#### Mode A: `realign_unwarp` (GRE fieldmap)

1. Decompress all `.nii.gz` to `.nii` (SPM cannot read compressed NIfTI)
2. For each task: calculate a VDM from the GRE fieldmap, matched to that task's first EPI volume
3. Realign & Unwarp all tasks in one multi-session SPM batch (joint motion + distortion correction)
4. Average per-task mean images into a session mean
5. Coregister T2w to session mean EPI; reslice ROI masks (nearest-neighbor interpolation)

**Output in func/:** `u<BOLD>.nii`, `meanu<BOLD>.nii` (one per task), `meanu_session.nii`, `rp_<BOLD>.txt`
**Output in anat/:** `r<ROI>.nii` (resliced masks)

#### Mode B: `topup` + `applytopup`

1. Decompress all `.nii.gz` to `.nii`
2. Extract first volume from the first task's BOLD (forward b0) and from the reverse-PE EPI (reverse b0)
3. Merge forward + reverse into a single 2-volume 4D file (`topup_merged_b0.nii.gz`)
4. Write `topup_acqparams.txt` with PE direction vectors and total readout time
5. Run FSL `topup` on the merged b0 to estimate the distortion field
6. For each task: run `applytopup --method=jac` to apply distortion correction to the full 4D BOLD
7. Realign (estimate + write) all corrected BOLDs in one SPM batch (motion correction only, since distortion is already corrected)
8. Coregister T2w to mean EPI; reslice ROI masks

**Output in func/:** `u<BOLD>_dc.nii`, `mean<BOLD>_dc.nii`, `rp_<BOLD>_dc.txt`, `topup_results_*`, `topup_acqparams.txt`
**Output in anat/:** `r<ROI>.nii`

#### Mode C: `topup` + `vdm`

1. Same as Mode B steps 1-5 (topup field estimation)
2. For each task: convert the topup Hz fieldmap to an SPM-compatible VDM using: `VDM = field_Hz * readout_sec * voxSize_PE`
3. Realign & Unwarp all tasks in one SPM batch using these VDMs (joint motion + distortion correction, same as `realign_unwarp` but with topup-derived VDMs)
4. Average per-task means into session mean; coregister + reslice ROIs

**Output in func/:** `u<BOLD>.nii`, `meanu<BOLD>.nii`, `meanu_session.nii`, `rp_<BOLD>.txt`, `vdm_task-*.nii`, `topup_results_*`
**Output in anat/:** `r<ROI>.nii`

#### Mode D: `realign_only`

1. Decompress all `.nii.gz` to `.nii`
2. Realign (estimate + write) all tasks in one SPM batch (motion correction only)
3. Coregister T2w to mean EPI; reslice ROI masks

**Output in func/:** `u<BOLD>.nii`, `mean<BOLD>.nii`, `rp_<BOLD>.txt`
**Output in anat/:** `r<ROI>.nii`

---

#### Motion Regressors (rp files)

All modes produce `rp_*.txt` files with 6 columns: **3 translations (mm) and 3 rotations (radians)**, one row per volume. What differs is what the motion parameters represent:

| Mode | rp file attached to | Meaning |
|------|---------------------|---------|
| `realign_unwarp` | `rp_<BOLD>.txt` | Motion estimated jointly with distortion correction |
| `topup` + `applytopup` | `rp_<BOLD>_dc.txt` | Motion estimated on already distortion-corrected images |
| `topup` + `vdm` | `rp_<BOLD>.txt` | Motion estimated jointly with VDM-based distortion correction |
| `realign_only` | `rp_<BOLD>.txt` | Pure rigid-body motion, no distortion model |

---

#### Preprocessing Mode Choice

- **realign_unwarp:** Uses GRE fieldmap + SPM Realign & Unwarp. Joint motion-distortion model. Requires correct readout time, TE1, TE2, and blip direction.

- **topup (applytopup):** Uses reverse-PE EPI + FSL topup. Distortion correction applied first (via `applytopup --method=jac`), then SPM Realign for motion. Requires correct PE directions and total readout time. Good when GRE fieldmaps are unavailable but reverse-PE EPIs exist.

- **topup (vdm):** Uses reverse-PE EPI + FSL topup, but converts the field to an SPM VDM so that SPM Realign & Unwarp handles both motion and distortion jointly. Same requirements as topup+applytopup.

- **realign_only:** Motion correction only. No fieldmap needed. Use if no fieldmap/reverse-PE is available.

**Input (all modes):**
- Native 4D fMRI (`sub-XX_ses-YY_task-X_bold.nii.gz`)
- Native T2w anatomical
- Native T2w ROI masks (left/right)
- Fieldmap files (mode-dependent, see above)

**SLURM Parameters (from config):**
- `PREPROC_CPUS` — CPUs per preprocessing job (default: 10)
- `PREPROC_MEM_PER_CPU` — Memory per CPU (default: 12G)
- `PREPROC_TIME` — Wall-clock limit per job (default: 10:00:00)
- `SLURM_PARTITION` — HPC partition name (default: compute)

---

#### Smoothing

The last step is an optional smoothing kernel you can set in the cfg file.
- Set to 0 if you want no smoothing
- Keep 5 mm FWHM if unsure or set to a custom smoothing kernel


**QC Checkpoint: Coregistration**

After preprocessing finishes, spot-check a few subjects:

1. Load unwarped EPI (`uxxxx.nii.gz`) and resliced ROI mask (`wroi_*.nii.gz`) in ITK-SNAP
2. Verify:
   - **Position:** Mask sits on expected anatomical signal (entorhinal cortex, hippocampus, etc.), not floating or shifted
   - **Coverage:** ROI not missing entire subregion due to coregistration error
   - **Distortion:** EPI anatomy not wildly stretched or banana-shaped

**Common Issues:**

| Issue | Likely Cause | Fix |
|-------|--------------|-----|
| EPI globally warped | Unwarping failed; wrong params (readout time, TE, deltaTE) | Re-check JSON parameters; try realign_only as baseline |
| MTL/EC signal absent | Over-aggressive unwarp, wrong fieldmap, wrong PE direction | Verify fmap corresponds to this run/session; inspect raw signal |
| Mask position off | Coregistration error, orientation mismatch, native T2w not native | Check T2w is truly native; verify coregistration reference |
| Multiple rp files per run | Previous preprocessing attempt lingering | Delete old derivatives; restart from native |

---

### Step 3: Prepare for GridCAT

```bash
./run_prep.sh
```

**What it does:**

Runs `prepare_gridcat_directory.m` for each subject-session:

1. **Split 4D → 3D:** Converts `uxxxx.nii.gz` into individual 3D volumes (GridCAT expects separate files per timepoint)
2. **Copy motion regressors:** Copies `rp_*.txt` to GLM directory with consistent naming
3. **Create bilateral ROIs:** Combines left and right masks into a bilateral mask (if `ROI_MODE=bilateral` in config)
4. **Verify directory structure:** Ensures GLM working directory exists and is populated

**Input:**
- Unwarped 4D EPI (`uxxxx.nii.gz`)
- Resliced ROI masks (`wroi_left.nii.gz`, `wroi_right.nii.gz`)
- Motion regressors (`rp_*.txt`)

**Output:**
- 3D EPI volumes: `GLM_runauto/sub-XX_ses-YY_task-X_run-Y_bold_0001.nii.gz`, `0002.nii.gz`, etc.
- Motion files: `GLM_runauto/rp_sub-XX_ses-YY_task-X_run-Y.txt`
- ROI masks: `GLM_runauto/sub-XX_ses-YY_roi_left.nii.gz`, `roi_right.nii.gz`, `roi_bilateral.nii.gz`

**Important:** Delete or isolate original 4D files after splitting. GridCAT must not accidentally mix 4D and 3D.

**QC Checkpoint: Bilateral ROIs**

Verify that `roi_bilateral.nii.gz` is the union of left and right without overlap (unless expected at midline).

---

### Step 4: Run GridCAT Analysis

```bash
./run_gridcat.sh
```

**What it does:**

Runs `run_gridcat_analysis.m`, which executes:

1. **GLM1:** Fits a single GLM per session using all requested tasks
   - Stacks all 3D volumes, motion regressors, and event timing
   - Estimates BOLD signal modulation by stimulus

2. **GLM2:** Runs directional connectivity analysis for each ROI
   - Extracts time series from each ROI
   - Computes voxel-wise connectivity (grid analysis) within the ROI
   - Outputs grid metrics (direction selectivity, hexadirectional modulation, etc.)

**Config Options (from pipeline_config.cfg):**

- **TASKS** — Which task numbers to include (e.g., `1,2,3`)
- **TR** — Repetition time (seconds)
- **HPF_CUTOFF** — High-pass filter cutoff (seconds), default 128
- **MASKING_THRESHOLD** — SPM masking threshold (0–1), higher = stricter signal mask
- **X_FOLD_SYMMETRY** — Grid symmetry (6 = hexagonal, 4 = square)
- **ROI_MODE** — Analyze `both`, `bilat_only`, or `lr_only`
- **MICROTIME_ONSET / MICROTIME_RESOLUTION** — SPM microtime settings

**Input:**
- 3D EPI volumes (from Step 3)
- Motion regressors
- EventData.txt files (task timing)
- ROI masks

**Output:**
- GridCAT results in `GLM_runauto/` and subdirectories
- Grid metrics per ROI:
  - Direction maps
  - Hexadirectional modulation maps
  - Rayleigh statistics
  - NaN/coverage maps
- Intermediate files (design matrices, activation maps, etc.)

**EventData.txt Preparation**

GridCAT requires event files with precise timing. Use an automated script (e.g., Poppy's converter) or create manually:

**Format:**
```
onset duration 1 condition_label
onset duration 1 condition_label
...
```

Where:
- `onset` = time in seconds from run start
- `duration` = stimulus/event duration in seconds
- `1` = event weight (usually 1)
- `condition_label` = text label for condition

**QC:**
- Verify event counts per condition match raw task logs
- Check that all runs are represented
- Confirm no temporal overlaps (unless intended)

---

### Step 5: Collect Results and Generate Report

```bash
Rscript collect_gridcat_output.R
```

**What it does:**

Reads all GridCAT output files and produces:

1. **CSV file** with grid metrics per subject-session-ROI
   - Direction selectivity index
   - Hexadirectional modulation strength
   - NaN/dropout percentage
   - Rayleigh uniformity test
   - Other relevant statistics

2. **HTML report** with:
   - Summary tables (overview of all subjects)
   - Diagnostic plots (QC metrics, coverage, voxel stability)
   - Subject-level detail pages

**Input:**
- GridCAT results (from Step 4)
- Pipeline config (for study metadata)

**Output:**
- `gridcat_results.csv` (for statistical analysis)
- `gridcat_report.html` (interactive overview)

**QC Interpretation:**

| Metric | Interpretation | Action if Problem |
|--------|-----------------|-------------------|
| High NaN % (>50%) | Poor signal or dropout in ROI | Check mask overlay; verify coregistration; inspect raw EPI signal |
| Low connectivity strength | Weak or absent hexadirectional modulation | Could be real (no grid signal) or pipeline artifact; verify event timing and motion regressors |
| Rayleigh p < 0.05 | Directional signal detected | Normal; suggests grid coding present |
| Rayleigh p > 0.05 | Uniform direction distribution | Could be real noise or no grid signal; check preprocessing quality |

---

## Using the Master Pipeline Runner

Instead of running steps individually, use the interactive menu:

```bash
./run_pipeline.sh
```

**Menu Options:**

```
GridCAT Pipeline - Main Menu
============================
0) Discover subjects/sessions (step0_make_subses_list.sh)
F) Filter subject/session list (filter_subses_list.sh)
1) Validate configuration (validate_config.sh)
2) Run SPM preprocessing (sbatch spm_preproc_array.sbatch)
3) Prepare for GridCAT (run_prep.sh)
4) Run GridCAT analysis (run_gridcat.sh)
5) Collect results (collect_gridcat_output.R)
A) Run all steps (0-5 in sequence)
Q) Quit

Select [0-5/F/A/Q]:
```

**Recommended workflow:**

1. Run `0` (discover) once at setup
2. Run `F` (filter) if you need to select specific subjects/sessions
3. Run `1` (validate) before first full run and after config changes
4. Run `2` (preprocessing) via HPC; wait for completion
5. Run `3` → `4` → `5` in sequence
6. Or run `A` to automate all steps (filter runs automatically after step 0)

**DRY_RUN Mode:**

Set `DRY_RUN=true` in config to preview actions without executing them.

```bash
./run_pipeline.sh
# Select step 6
# Scripts will print intended actions but not execute
```

Useful for validating configuration and checking file paths before committing to long HPC runs.

---

## Mask Creation (T2w Space)

### What You Need

- Left and right ROI masks in native T2w space
- Clean, contiguous labels (integer-valued)
- Anatomically accurate segmentation

### Options

1. **Manual Segmentation** — Draw masks using ITK-SNAP or similar; follow standard anatomical protocol
2. **Automated Tools** — Use ASHS (Automatic Segmentation of Hippocampal Subfields) if available
   - Use `ashs_extract_config.sh` and `ashs_extract_labels.sh` helpers to extract ROI outputs

### QC Checkpoint: Mask Anatomy

1. Open T2w and ROI masks in ITK-SNAP
2. Verify:
   - **Anatomical plausibility:** Mask follows expected anatomical boundaries
   - **No left/right swap:** Check orientation and compare to known templates
   - **No excess CSF/WM:** ROI should not include large ventricular or white matter chunks unless intended
   - **Label consistency:** If using multiple ROI subregions, verify labeling is consistent

### Typical Failure Modes

| Problem | Cause | Fix |
|---------|-------|-----|
| Mislabeled hemisphere | Orientation convention mismatch or manual error | Verify against template; swap if needed |
| Different orientations (T2w vs masks) | Masks created in different software/space | Reorient masks to match T2w using SPM or FSL |
| Non-binary or fuzzy masks | Probability maps instead of label masks | Use thresholding or rerun segmentation with label output |
| ROI outside brain | Manual error or space mismatch | Re-check in ITK-SNAP; redo if necessary |

---

## Preprocessing: Choosing a Mode

### Decision Tree

**Use `realign_only` if:**
- No fieldmap or reverse-PE EPI available
- Study has minimal geometric distortion
- You want simplicity and faster processing
- Diagnostic step to rule out unwarping issues

**Use `realign_unwarp` if:**
- GRE fieldmap (magnitude1 + phasediff) is available and high-quality
- You have verified acquisition parameters (readout time, TE1, TE2, blip direction)
- MTL/EC signal is critical and distortion is significant

**Use `topup` if:**
- Reverse phase-encode EPI is available (but no GRE fieldmap, or you prefer topup)
- Works with either `applytopup` or `vdm` sub-method (see below)

### Option A: `realign_only`

```
PREPROC_MODE=realign_only
```

- Runs SPM Realign (motion estimation only)
- Outputs: `rp*.txt` (6 motion parameters per timepoint), `u*.nii` (resliced)

**Pros:** Simple, no parameter tuning, faster
**Cons:** Distortion remains uncorrected, especially near skull base

### Option B: `realign_unwarp`

```
PREPROC_MODE=realign_unwarp
```

- Computes VDM from GRE fieldmap per task
- Runs SPM Realign & Unwarp (joint motion + distortion correction)
- Outputs: `u*.nii`, `meanu*.nii`, `rp*.txt`

**Pros:** Joint motion-distortion model; motion parameters account for distortion interactions
**Cons:** Parameter-sensitive; requires correct readout time, TE1, TE2

### Option C: `topup` (with sub-methods)

```
PREPROC_MODE=topup
TOPUP_APPLY_METHOD=applytopup   # or: vdm
```

Uses FSL topup to estimate the distortion field from a pair of opposite phase-encode EPIs.

**Sub-method `applytopup`:**
- FSL `applytopup --method=jac` corrects distortion first, then SPM Realign handles motion
- Output BOLDs have `_dc` suffix (e.g., `u<SUB>_<SES>_task-run1_bold_dc.nii`)
- Motion parameters are estimated on the already corrected images

**Sub-method `vdm`:**
- Topup Hz fieldmap is converted to an SPM VDM
- SPM Realign & Unwarp handles motion + distortion jointly (same as `realign_unwarp`, but with topup-derived VDMs)
- Output naming is the same as `realign_unwarp`

**Pros:** Works without GRE fieldmaps; topup is well-validated for EPI distortion correction
**Cons:** Requires a reverse-PE EPI; `applytopup` method separates motion and distortion into two steps

#### Critical Parameters for `topup`

| Config key | What it is | Where to find it |
|-----------|-----------|-----------------|
| `TOPUP_PE_DIR_BOLD` | Phase-encode direction of task EPIs (FSL convention: x, x-, y, y-, z, z-) | BIDS JSON `PhaseEncodingDirection`: i=x, i-=x-, j=y, j-=y- |
| `TOPUP_PE_DIR_REVERSE` | Phase-encode direction of reverse-PE EPI (opposite of BOLD) | Same JSON field on reverse-PE file |
| `TOPUP_READOUT_SEC` | Total readout time in seconds | BIDS JSON `EstimatedTotalReadoutTime` (already in seconds) |
| `TOPUP_REVERSE_PE_PATTERN` | Filename pattern for the reverse-PE EPI | e.g., `task-reverse_bold` or `dir-PA_epi` |

**IMPORTANT:** `TOPUP_READOUT_SEC` may differ from `TOTAL_READOUT_MS` (the latter is for GRE fieldmaps, in milliseconds). Double-check both values independently from the BIDS JSON sidecars.

#### Critical Parameters for `realign_unwarp`

| Config key | What it is | Where to find it |
|-----------|-----------|-----------------|
| `TOTAL_READOUT_MS` | Total EPI readout time in milliseconds | BIDS func JSON `TotalReadoutTime` (multiply by 1000) |
| `TE_SHORT_MS` / `TE_LONG_MS` | Echo times of the GRE fieldmap | BIDS fmap JSON `EchoTime1` / `EchoTime2` (multiply by 1000) |
| `BLIP_DIRECTION` | Blip direction for SPM FieldMap (+1 or -1) | BIDS JSON `PhaseEncodingDirection`: j = +1, j- = -1 |

**Verification (for either mode):**

```bash
# For topup: check EstimatedTotalReadoutTime
jq '.EstimatedTotalReadoutTime' sub-XX_ses-YY_task-1_bold.json

# For realign_unwarp: check TotalReadoutTime and echo times
jq '.TotalReadoutTime' sub-XX_ses-YY_task-1_bold.json
jq '.EchoTime1, .EchoTime2' sub-XX_ses-YY_phasediff.json

# Phase-encode direction (both modes)
jq '.PhaseEncodingDirection' sub-XX_ses-YY_task-1_bold.json
```

---

## HPC Workflow

### Submitting Preprocessing to SLURM

**Prerequisite:** `subses_list.txt` exists (from Step 0)

**Method 1: Via run_pipeline.sh**

```bash
./run_pipeline.sh
# Select option 2
# Script prompts for HPC parameters, submits sbatch
```

**Method 2: Direct sbatch**

```bash
N=$(wc -l < subses_list.txt)
sbatch --array=1-$N spm_preproc_array.sbatch
```

### Monitoring Jobs

```bash
# Check queue
squeue -u $USER

# Check specific job
squeue -j JOBID

# Check completed job (last 100 lines)
tail -100 slurm-JOBID.out
```

### Common SLURM Issues

| Issue | Solution |
|-------|----------|
| "QOSMaxCpuPerUserLimit" | Reduce `PREPROC_CPUS` or wait for cluster availability |
| Job times out | Increase `PREPROC_TIME` in config; profile MATLAB runtime |
| MATLAB crashes | Increase `PREPROC_MEM_PER_CPU`; check for disk space |
| Array job not starting | Verify `subses_list.txt` is not empty; check cluster load |

### Checking Job Logs

Each job creates a log: `slurm-JOBID_array-TASKID.out`

Inspect for:
- MATLAB errors or warnings
- File not found messages
- Path issues
- Memory/time limit reached

### After Preprocessing

Once all jobs complete:

1. Spot-check a few subject outputs for coregistration quality (see QC Checkpoint above)
2. Verify `rp_*.txt` files exist for all runs
3. Proceed to Step 3 (Prepare for GridCAT)

---

## Coregistration: Bringing Masks into EPI Space

### Process

1. **Coregister T2w → EPI mean** (not the reverse; preserves functional space)
2. **Apply transform to ROI masks** using nearest-neighbor interpolation (label-safe)
3. **Optional:** Apply transform to T2w itself for QC overlays

### Critical Rule

**Always coregister T2w to EPI, not EPI to T2w**, to avoid distorting functional space.

### Interpolation

- **ROI masks:** Use nearest-neighbor (`interp=0`) or label-safe mode (`interp=-1` in SPM)
  - Preserves integer labels
  - Avoids smoothing mask boundaries
  - Output will look "blocky" (normal)

- **T2w reference:** Linear interpolation acceptable
  - Used only for visualization/QC
  - Not used in downstream analysis

### QC Checkpoint: Mask-EPI Overlay

After Step 2 preprocessing completes, verify:

1. Load `meanU*.nii.gz` (unwarped EPI mean) and `wroi_*.nii.gz` (resliced masks) in ITK-SNAP
2. Check:
   - **Position:** Mask overlaps expected anatomy (entorhinal cortex, hippocampus, etc.)
   - **No gross misalignment:** Mask not shifted by >5 mm
   - **Coverage:** Entire ROI visible in EPI (no clipping due to distortion)
   - **No anatomy mismatch:** Mask follows anatomy, not floating in CSF/white matter

3. Typical views:
   - Axial: Verify left/right orientation, mask within MTL
   - Coronal: Verify anterior-posterior position along long axis
   - Sagittal: Verify mask follows expected medial-lateral anatomy

### Troubleshooting Coregistration

| Symptom | Likely Cause | Action |
|---------|--------------|--------|
| Mask position off by >5 mm | Coregistration failure | Verify EPI reference (should be unwarped mean); check if T2w is truly native |
| Mask outside brain | Space mismatch or wrong reference | Check coregistration reference file; verify T2w native space |
| Mask rotated or sheared | Orientation/header issue | Check sform/qform consistency; reorient T2w if needed |
| Multiple rp files per run | Previous preprocessing lingering | Delete old derivatives; restart from native EPI |

---

## Preparing Data for GridCAT

### Event Data (Task Timing)

GridCAT requires precise event timing in `EventData.txt` files.

**Creation:**
- Use automated script (e.g., Poppy's converter) to generate from raw task logs
- Or manually create if raw logs are well-structured

**Format (plain text):**
```
5.2 1.0 1 face
6.5 1.0 1 house
8.0 1.0 1 face
...
```

Columns: `onset_sec duration weight condition_label`

**File Naming:**
- Must include sub, ses, task identifiers matching pipeline naming
- Must be placed in `GLM_runauto/` directory
- One file per session (or per task, depending on GLM strategy)

**QC:**
- Verify event counts match raw task logs
- Check for temporal overlaps (unless multi-event trials intended)
- Ensure all tasks/runs are represented
- Compare timing across subjects/runs for consistency

### Motion Regressors

Step 3 (Prepare) automatically copies motion regressors from preprocessing. Verify:

- **One rp per run:** No duplicates or missing files
- **Correct association:** Motion file corresponds to the EPI it was computed from
- **Correct location:** Copied to `GLM_runauto/` with consistent naming

**If mismatch suspected:**

```bash
# List all rp files
find derivatives/ -name "rp_*" -exec ls -l {} \;

# Check timestamps; newest should correspond to current preprocessing
# Delete old rp files if multiple exist per run
```

### 4D → 3D Splitting

Step 3 (Prepare) calls `prepare_gridcat_directory.m`, which:

1. Loads each 4D EPI file
2. Splits into individual 3D volumes
3. Writes: `sub-XX_ses-YY_task-X_run-Y_bold_NNNN.nii.gz` (NNNN = 0001, 0002, ...)

**Important:** After splitting, delete or isolate original 4D files:

```bash
# Move 4D files out of GLM directory to avoid mixing
mkdir GLM_runauto/4D_backups
mv GLM_runauto/*_bold.nii.gz GLM_runauto/4D_backups/
```

GridCAT expects 3D inputs only.

### Bilateral ROIs

If using bilateral analysis (`ROI_MODE=bilateral` in config):

- Step 3 automatically creates `roi_bilateral.nii.gz` as union of left and right
- Verify in ITK-SNAP that bilateral = left OR right (no gaps, minimal overlap)

---

## Running GridCAT GLM Analysis

### GLM1: Session-Level Model

- Fits single design matrix per session
- Includes all task runs with their respective event timings
- Models BOLD signal modulation by stimulus

### GLM2: ROI-Level Connectivity

- Extracts mean time series from ROI
- Computes voxel-wise correlation/directional connectivity within ROI
- Outputs grid metrics (direction selectivity, hexadirectional modulation, etc.)

### Key Settings

| Setting | Options | Notes |
|---------|---------|-------|
| **TASKS** | Comma-separated task numbers | E.g., `1,2,3`; exclude others via EXCLUDE_TASKS |
| **TR** | Seconds, e.g., 2 | Must match acquisition |
| **HPF_CUTOFF** | Cutoff in seconds | Default 128; set 0 to disable |
| **MASKING_THRESHOLD** | 0–1 | Higher = stricter signal mask; default 0.8 |
| **X_FOLD_SYMMETRY** | 6 (hexagonal) or 4 (square) | 6 for grid cells |
| **ROI_MODE** | `both`, `bilat_only`, or `lr_only` | Analyze hemispheres separately, combined, or both |
| **DERIVATIVES** | `0,0` / `1,0` / `1,1` | HRF derivatives: none / temporal / temporal+dispersion |

### Typical Configuration

```cfg
TASKS=1,2,3
TR=2
HPF_CUTOFF=128
MASKING_THRESHOLD=0.8
X_FOLD_SYMMETRY=6
ROI_MODE=both
DERIVATIVES=0,0
MICROTIME_ONSET=8
MICROTIME_RESOLUTION=16
```

### Output Structure

```
GLM_runauto/
  sub-XX_ses-YY_grid_metrics.mat         (MATLAB results)
  sub-XX_ses-YY_direction_map.nii.gz     (voxel-wise direction map)
  sub-XX_ses-YY_hexmod_map.nii.gz        (hexadirectional modulation map)
  sub-XX_ses-YY_rayleigh_stats.txt       (Rayleigh uniformity test)
  sub-XX_ses-YY_coverage.nii.gz          (voxel inclusion map)
```

---

## Interpreting Results

### Output Files

After Step 5 (collect_gridcat_output.R):

- **gridcat_results.csv** — Summary metrics per subject-session-ROI (for statistics)
- **gridcat_report.html** — Interactive visualization of QC metrics and results

### Key Metrics

| Metric | Interpretation |
|--------|-----------------|
| **Direction Selectivity** | Voxel's preference for movement direction (0-1, higher = more selective) |
| **Hexadirectional Modulation** | Strength of 6-fold directional symmetry (grid cell-like property) |
| **NaN %** | Fraction of voxels excluded due to low signal/thresholds (aim <30%) |
| **Coverage** | % of ROI voxels included in analysis (aim >70%) |
| **Rayleigh p-value** | Statistical test for uniformity of direction distribution (p<0.05 = directional signal) |
| **Voxel Stability** | Consistency of metrics across repetitions/tasks (if applicable) |

### Interpreting High NaN %

High NaN fractions (>50%) indicate problems:

**Check:**
1. **Mask overlay quality** — Verify mask properly coregistered to EPI (see QC Checkpoint, Coregistration section)
2. **Raw EPI signal** — Inspect original `uxxxx.nii.gz` in MTL region for signal dropout or artifact
3. **Thresholds** — Review `GRID_THRESHOLD` and other cutoffs in config; may be too strict
4. **Space/reslicing** — Confirm masks resliced with correct interpolation (nearest-neighbor for labels)

**Fix:**
- Re-run coregistration if alignment is off
- Relax thresholds in config if EPI quality is acceptable
- Exclude problematic runs or subjects from downstream analysis

### Interpreting Weak Connectivity

Weak or absent hexadirectional modulation could be:

1. **Real:** No grid signal in this ROI/subject/task (valid finding)
2. **Artifact:**
   - Wrong event timing in EventData.txt (timing off, events missing)
   - Wrong motion regressors (using rp from different preprocessing)
   - Poor signal-to-noise ratio due to distortion, dropout, or artifact
   - Overly strict GLM thresholds

**Verify:**
- Event counts per condition match raw task logs
- Motion regressors correspond to final (unwarped) EPI
- Check raw EPI signal quality (no dropout, reasonable SNR)
- Inspect design matrix (events properly modeled)

### Rayleigh Test Interpretation

- **p < 0.05:** Evidence of directional modulation (grid-like signal detected)
- **p > 0.05:** No significant directional preference (either no grid signal or uniform noise)

Both outcomes are valid; interpretation depends on your research question.

---

## Troubleshooting

### Common Issues and Solutions

#### Issue: "File not found" during preprocessing

**Cause:** BIDS naming mismatch, missing file, or path error

**Fix:**
1. Verify dataset path in `pipeline_config.cfg`
2. Check BIDS structure (especially anat, fmap, func directories)
3. Run `validate_config.sh` to diagnose
4. Inspect `subses_list.txt` for correct subject/session format
5. Use `ls` and `find` to locate missing files manually

#### Issue: "MATLAB error" or script crashes

**Cause:** Insufficient memory, MATLAB path issues, or license problem

**Fix:**
1. Increase `HPC_MEMORY` in config (32→64 GB)
2. Check MATLAB license availability: `matlab -r "license" -nodisplay`
3. Verify SPM path in `read_pipeline_config.m`
4. Check disk space: `df -h`
5. Check slurm log for specific error message

#### Issue: Coregistration looks wrong

**Cause:** Bad reference image, orientation mismatch, or native T2w preprocessed

**Fix:**
1. Verify EPI reference is unwarped mean (`meanU*.nii.gz`), not native
2. Check T2w is truly native (not previously resliced/coregistered)
3. Inspect sform/qform consistency: `nifti_tool -disp_hdr file.nii.gz | grep -A5 sform`
4. Check orientation in ITK-SNAP (verify left/right labels)
5. Re-run coregistration with explicit reference image if needed

#### Issue: Unwarping/topup looks wrong (EPI distorted, signal missing)

**Cause:** Incorrect readout time, TE, DeltaTE, or PE direction parameters

**Fix:**
1. Re-derive total readout time from BIDS JSON: check `TotalReadoutTime` (for GRE) or `EstimatedTotalReadoutTime` (for topup)
2. For `realign_unwarp`: verify TE1, TE2 from fmap JSON
3. For `topup`: verify `TOPUP_PE_DIR_BOLD` and `TOPUP_PE_DIR_REVERSE` match the BIDS `PhaseEncodingDirection` field (i=x, j=y, etc.)
4. Confirm PE direction is correct (swapped PE = correction applied in wrong direction)
5. Cross-check with original acquisition protocol (ask MRI tech if unsure)
6. Try `realign_only` as diagnostic baseline to confirm preprocessing is otherwise correct
7. Check `topup_acqparams.txt` in the func/ directory — the two lines should have opposite PE vectors and the readout time should match your BIDS JSON

#### Issue: Multiple rp files per run

**Cause:** Previous preprocessing attempt not cleaned up

**Fix:**
1. Identify which rp files are current: `ls -l rp_*.txt | sort -k6,7`
2. Newest timestamps = current preprocessing; delete others
3. Or delete entire derivatives/ folder and restart from native EPI

#### Issue: GridCAT output is empty or all NaN

**Cause:** Mask not in correct space, missing event data, or GLM failed

**Fix:**
1. Verify 3D EPI files exist in `GLM_runauto/` (from Step 3)
2. Verify EventData.txt has correct format and timing
3. Verify ROI masks are in `GLM_runauto/` with correct naming
4. Check GLM log for specific error message
5. Verify motion regressors (rp*.txt) are copied and accessible

#### Issue: High NaN % in results (>50%)

**Cause:** Poor coregistration, signal dropout, or overly strict thresholds

**Fix:**
1. Re-check mask-EPI overlay in ITK-SNAP (see QC Checkpoint, Coregistration)
2. Inspect raw EPI for signal dropout or artifact in MTL region
3. Relax `GRID_THRESHOLD` in config (try 0.10 instead of 0.05)
4. Verify mask resliced with nearest-neighbor (not linear) interpolation
5. Check mask is not outside brain or rotated

---

## Study Design Decisions

These decisions must be made explicitly for each project and documented in the config:

### Task Inclusion/Exclusion

- Which tasks are GridCAT-relevant?
- Exclude resting-state or other conditions?
- How to handle multiple runs of same task (if present)?

### Preprocessing Mode

- `realign_unwarp` (GRE fieldmap), `topup` (reverse-PE EPI), or `realign_only`?
- If using `topup`: `applytopup` (two-step) or `vdm` (joint model)?
- Justified by acquisition quality, fieldmap availability, and prior experience

### ROI Definition

- Which anatomical structures (hippocampus, entorhinal cortex, etc.)?
- Bilateral combined or hemispheres separate?
- Automated (ASHS) or manual segmentation?

### GLM Design

- Which regressors (motion, CSF, white matter)?
- Smoothing kernel size? High-pass filter?
- Event model (boxcar, stick, parametric)?

### QC Thresholds

- NaN tolerance (e.g., <30%)
- Voxel stability requirement
- Dropout criteria (e.g., exclude if >50% NaN)

All of these are configurable in `pipeline_config.cfg` and should be documented in your project's methods section.

---

## Key Principles

1. **Single config file controls everything.** Users edit `pipeline_config.cfg` once per project; all scripts read from it automatically.

2. **Master runner for convenience.** Use `run_pipeline.sh` to execute steps via interactive menu, or run steps individually.

3. **Auto-detection of paths.** Each script determines its own directory and finds `pipeline_config.cfg` automatically; no manual path editing needed.

4. **Full user control.** All settings are visible and adjustable; no hidden defaults or black-box processing.

5. **Pipeline works on any BIDS dataset.** Change only the config to process a new study.

6. **Each step can be run independently.** No hard dependency on step order, though recommended sequence is 0→1→2→3→4→5.

7. **DRY_RUN mode for preview.** Set `DRY_RUN=true` to see intended actions without executing.

8. **Fail gracefully on missing data.** Set `FAIL_ON_MISSING=false` for lenient processing; set to `true` to halt on any missing file.

---

## Final Checklist Before Analysis

- [ ] BIDS dataset validated (all required files present)
- [ ] `pipeline_config.cfg` edited with project-specific settings
- [ ] `validate_config.sh` passes all checks
- [ ] `subses_list.txt` contains all expected subject-sessions
- [ ] Subject/session filters configured (if needed) and `filter_subses_list.sh` run
- [ ] Fieldmap/topup parameters verified against BIDS JSONs (readout time, TE, PE direction)
- [ ] ROI masks created and QC'd in T2w space
- [ ] EventData.txt files created for all tasks
- [ ] SPM preprocessing submitted to HPC and monitoring logs
- [ ] Coregistration quality verified on sample subjects
- [ ] 4D → 3D splitting completed; old 4D files removed from GLM directory
- [ ] GridCAT GLM analysis completed
- [ ] Results collected and report generated
- [ ] Output CSV and HTML reviewed for anomalies

---

## Support & Questions

If issues arise:

1. Check the Troubleshooting section above
2. Verify config settings match your study design
3. Inspect logs: `slurm-*.out` for HPC jobs, MATLAB command window output for script errors
4. Verify BIDS naming and file locations manually
5. Run `validate_config.sh` to diagnose configuration problems

Document any modifications to scripts or config; keep a record of preprocessing decisions for reproducibility and future reference.

---

**End of SOP** 
