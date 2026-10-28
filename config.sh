# =============================================================================
#  Grid Cell Pipeline — shell settings
# =============================================================================
#
#  Paths, cluster resources and which subjects/runs to process.
#  Sourced directly by the shell scripts, so it is plain bash:
#  KEY=value, no spaces around '=', no quotes, no trailing spaces.
#
#  SPM and GridCAT parameters live in spm_options.m, not here.
#
#  To see your current settings:   cat config.sh
# =============================================================================


# -----------------------------------------------------------------------------
#  Paths
# -----------------------------------------------------------------------------

# BIDS data. Read only — the pipeline never writes here.
BIDS_ROOT=/sc-projects/sc-proj-cc02-brainspace/myspace/julius/SOLIS_GC/bids

# Where output goes. Created if missing.
OUTPUT_ROOT=/sc-projects/sc-proj-cc02-brainspace/myspace/julius/SOLIS_GC/workDir

# This script directory on the cluster.
PIPELINE_DIR=/sc-projects/sc-proj-cc02-brainspace/myspace/julius/SOLIS_GC/Grid-Cell-Pipeline

SPM_DIR=/sc-projects/sc-proj-cc02-brainspace/software/matlab/spm25
GRIDCAT_DIR=/sc-projects/sc-proj-cc02-brainspace/myspace/julius/sirius_wip/GridCAT

# Ships with GridCAT. Empty = GRIDCAT_DIR/CircStat2012a
CIRCSTAT_DIR=


# -----------------------------------------------------------------------------
#  Derivatives — everything the pipeline writes
# -----------------------------------------------------------------------------
#
#  SPM's Realign & Unwarp and Coregister write the transforms they estimate
#  into the headers of their INPUT images, and reading a .nii.gz leaves a
#  decompressed twin behind. So every input is copied here first and the copy
#  is what gets touched. Run any of that straight on BIDS_ROOT and the raw
#  data quietly stops being raw.

# Empty = <OUTPUT_ROOT>/derivatives
DERIV_ROOT=

# Three datasets inside it, each a BIDS derivatives dataset of its own.
#   spm-preproc  regenerated from scratch on every full run (see DERIV_RESET)
#   rois         kept: manual segmentations, re-importing them is your call
#   fieldmaps    kept: rebuilding means another FSL run
DERIV_PREPROC=spm-preproc
DERIV_ROIS=rois
DERIV_FIELDMAPS=fieldmaps

# When to clear the preprocessing derivatives:
#   auto    clear when the run starts preprocessing from its beginning.
#           A resumed run keeps what is there, because that is its input.
#   always  clear on every run. Resuming stops working.
#   never   never clear; new output overwrites old file by file.
DERIV_RESET=auto

# Names this run's GridCAT output folder: empty -> GLM_output,
# otherwise GLM_output_<variant>. Each submission snapshots the config, so
# you can edit this and submit the next variant while one is still running.
# Letters, digits, underscore, dash and dot only.
RUN_VARIANT=6f_t5_am_keep1_avg_1


# -----------------------------------------------------------------------------
#  What to process
# -----------------------------------------------------------------------------

# Task labels, comma-separated, matching task-<LABEL> in the BIDS filenames.
TASKS=1,2,3

# How to correct distortion:
#   realign_unwarp   GRE fieldmap: needs phasediff + magnitude1 in fmap/.
#                    SPM unwraps the phase itself.
#   topup            FSL topup from a reverse-PE EPI.
#   precalc_fieldmap ready-made B0 fieldmap, built by tools/make_fieldmaps.sh
#                    for sessions that have a phasediff but NO magnitude.
#   realign_only     motion correction only, no distortion correction.
PREPROC_MODE=realign_unwarp

# Which preprocessing stages to run, comma-separated. Each reads what it
# needs from disk, so you can stop after any of them, check or replace the
# intermediate files, and resume later.
#
#   topup     FSL topup distortion field.  -> func/topup_results_*
#             Does nothing unless PREPROC_MODE=topup.
#   vdm       one voxel displacement map per task.  -> func/vdm_task-<N>.nii
#             Does nothing for realign_only, or topup+applytopup.
#   slicetime shift every slice to the same point in the TR.  -> a*_bold.nii
#             Does nothing when DO_SLICE_TIMING=0.
#   realign   apply the correction + motion correction.
#             -> u*_bold.nii, rp_*.txt, meanu_session.nii
#   coreg     coregister T2w -> session mean, reslice ROIs.  -> anat/r*_mask.nii
#   smooth    Gaussian smoothing.  -> func/su*_bold.nii
#
# Shorthands: all, fieldmap (= topup,vdm), post_fieldmap (= slicetime,realign,coreg,smooth)
# Override for one submission:  bash 01_preprocess.sbatch --stages smooth
PREPROC_STAGES=all

# Coregister T2w -> EPI and reslice ROI masks into EPI space? (1=yes, 0=skip)
DO_COREG_ROIS=1

# Coregister the T2w to the session mean THROUGH a whole-brain EPI, named
# here by its task label ('rest' finds func/*_task-rest*_bold.nii). Empty
# means the direct route.
#
# Worth turning on when the T2w and the task EPI are both partial slabs. They
# then overlap in only a small block of tissue, and normalised mutual
# information over a small overlap has a degenerate optimum: pushing the
# images apart shrinks the overlap to a near-uniform sliver that can score
# better than the true alignment. sub-05 ses-01 failed exactly this way -
# the T2w ended up 141 mm away and both ROI masks resliced to zero voxels.
# A whole-brain EPI overlaps both images generously, so each step is easy.
#
# The masks are still resampled only ONCE: both coregistration steps are
# estimate-only, and the single reslice at the end applies the composed
# transform. Two rigid transforms compose into one rigid transform, so this
# costs no accuracy over the direct route.
COREG_VIA_TASK=

# Fail the job when a resliced ROI mask comes back with fewer than this many
# voxels. A coregistration that diverges does not raise an error - it writes
# a mask that misses the brain, and the session then disappears from the
# analysis with nothing to show why. 0 disables the check.
COREG_MIN_ROI_VOXELS=5


# -----------------------------------------------------------------------------
#  Subject / session filtering
# -----------------------------------------------------------------------------
#  Applied by 00_list_subjects.sh. All optional; empty = include everything.
#  INCLUDE is a whitelist and is applied first, then EXCLUDE removes from it.
#
#   INCLUDE_SUBJECTS=sub-01s13,sub-02s13
#   EXCLUDE_SESSIONS=ses-03,ses-05
INCLUDE_SUBJECTS=
INCLUDE_SESSIONS=
EXCLUDE_SUBJECTS=sub-07s17
EXCLUDE_SESSIONS=

# Skip subject-sessions that have no ROI masks in anat/ (saves time).
REQUIRE_ROI=true

# Tasks to exclude from the GridCAT analysis, comma-separated.
EXCLUDE_TASKS=rest,reverse,mean


# -----------------------------------------------------------------------------
#  ROI masks
# -----------------------------------------------------------------------------

# Filenames the pipeline looks for in anat/. Extra BIDS entities tolerated.
ROI_PATTERN_LEFT=_hemi-left_label-lpmErC_mask.nii
ROI_PATTERN_RIGHT=_hemi-right_label-rpmErC_mask.nii

# The anatomical label in those filenames. Everything downstream — the
# GridCAT prep and the GLM — filters on this, so it must match the patterns.
ROI_LABEL=ErC

# Where tools/move_rois.sh imports manual segmentations from. Files are
# expected as: sub-XX_ses-YY_[run-N_]<ROI_LABEL>_<hemi>.nii[.gz]
ROI_SOURCE_DIR=/sc-projects/sc-proj-cc02-brainspace/myspace/julius/BRACE_GC_testing/v2_pipe_t5_data/ROI_lucie/ErC_ROIs_BIDS

# ROI analysis mode for GridCAT GLM2: both | bilat_only | lr_only
ROI_MODE=both


# -----------------------------------------------------------------------------
#  Multiple runs of the same scan
# -----------------------------------------------------------------------------
#
#  Run numbers are counted per modality and say nothing across modalities:
#  if the subject climbed out of the scanner between sequences, run-2 of the
#  T2w and run-2 of a task are from different visits, and pairing them
#  coregisters the wrong anatomy. So the pipeline decides in this order:
#    1. your choice in RUN_SELECTION_FILE
#    2. the run acquired closest in time, from AcquisitionTime in the JSONs
#    3. nothing else is safe — see STRICT_RUN_MATCHING

# Written by: bash tools/scan_multirun.sh
# Empty = rely on acquisition-time matching alone.
RUN_SELECTION_FILE=/sc-projects/sc-proj-cc02-brainspace/myspace/julius/SOLIS_GC/Grid-Cell-Pipeline/run_selection.tsv

# How far apart two scans may be acquired and still count as the same visit
# to the scanner, in minutes.
RUN_MATCH_GAP_MIN=20

# When nothing can decide between several runs:
#   true   stop and say so. Silently pairing the wrong T2w with a task is a
#          result you cannot spot afterwards, so this is the default.
#   false  take the highest run number and warn loudly.
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

# Filename substring for the reverse-PE EPI. Examples: dir-RL_epi,
# dir-PA_epi, task-reverse_bold
TOPUP_REVERSE_PE_PATTERN=task-reverse_bold

# FSL topup config file. Empty = FSL default b02b0.cnf
TOPUP_CONFIG=

# How to apply the correction:
#   vdm         convert the topup field to an SPM VDM, then Realign & Unwarp
#   applytopup  correct the EPIs with FSL applytopup, then SPM Realign only
TOPUP_APPLY_METHOD=vdm

# Pre-computed topup output, used when the "topup" stage is NOT run. Give
# the prefix without the _fieldcoef.nii.gz / _movpar.txt / _field.nii.gz
# suffix. {SUB} and {SES} are substituted. Empty = use func/topup_results.
TOPUP_EXISTING_PREFIX=

# --- precalc_fieldmap ---
# Patterns for fmap/, WITHOUT the extension (unlike ROI_PATTERN_* above).
FIELDMAP_PATTERN=_fieldmap
FIELDMAP_MAGNITUDE_PATTERN=_magnitude

# Units on disk. fsl_prepare_fieldmap writes rad/s; SPM wants Hz, so rad/s
# is divided by 2*pi before use. Getting this wrong scales the whole
# distortion correction by 2*pi.
#   rad/s = output of fsl_prepare_fieldmap (default)
#   Hz    = already converted
FIELDMAP_UNITS=rad/s


# -----------------------------------------------------------------------------
#  Building fieldmaps from a structural (tools/make_fieldmaps.sh)
# -----------------------------------------------------------------------------
#  For sessions with a phase-difference image but no magnitude, so
#  fsl_prepare_fieldmap has nothing to unwrap against. Builds the missing
#  magnitude from the T1w (bet -> flirt into the phasediff grid -> mask).
#  Run it with --dry-run first, then set PREPROC_MODE=precalc_fieldmap.

# Passed straight to fsl_prepare_fieldmap: SIEMENS or PHILIPS
FIELDMAP_SCANNER=SIEMENS

# Echo time difference in ms. Empty = read per session from the phasediff
# JSON (EchoTime2 - EchoTime1), which is safest.
FIELDMAP_DELTA_TE=

# bet fractional intensity threshold (lower keeps more brain).
FIELDMAP_BET_F=0.5

# Where the magnitude comes from:
#   structural  always build it from the T1w. Every session is processed
#               identically, which matters when you compare sessions.
#   auto        use the real *_magnitude1 when there is one, T1w otherwise.
#               Better per session, mixes two kinds across a longitudinal set.
FIELDMAP_MAGNITUDE_SOURCE=structural


# -----------------------------------------------------------------------------
#  Output file prefixes
# -----------------------------------------------------------------------------
#  These must match what SPM actually writes. FUNC_PREFIX is derived from
#  them (SMOOTH_PREFIX + RESLICE_PREFIX when smoothing is on), so there is
#  nothing to keep in sync by hand.

RESLICE_PREFIX=u
SMOOTH_PREFIX=s
ROI_PREFIX=r
SLICETIME_PREFIX=a

# Correct for the fact that slices within a volume are acquired at
# different times? (1=yes, 0=no)
#
# Needs SliceTiming in the BOLD JSON sidecars - the stage reads the actual
# times per run rather than assuming an order, because correcting with the
# wrong order is worse than not correcting. It fails loudly if SliceTiming
# is absent instead of guessing.
#
# Worth it when slices span a large part of the TR. Here 30 slices are
# acquired sequentially over 1.97 s of a 2.06 s TR, so two slices can be
# almost a whole TR apart.
#
# Changing this changes which files every later stage reads (the corrected
# ones carry SLICETIME_PREFIX), so re-run the preprocessing from "slicetime"
# onwards after changing it.
DO_SLICE_TIMING=1

# Smoothing kernel FWHM in mm, isotropic. 0 skips smoothing entirely.
#
# This lives here rather than in spm_options.m because it decides which
# files every later stage reads: with smoothing on, the GridCAT prep picks
# up su*_bold.nii, with it off, u*_bold.nii. That is a pipeline switch, not
# just an SPM parameter.
SMOOTH_FWHM=5


# -----------------------------------------------------------------------------
#  SLURM
# -----------------------------------------------------------------------------

# Check your HPC with: sinfo -s
SLURM_PARTITION=compute

# Job notification email. Empty disables it. When set you get two mails per
# full run_all.sh chain: one when all preprocessing array tasks finish, one
# when GridCAT finishes. FAIL notifications are sent too.
NOTIFY_EMAIL=julius.schaut@charite.de

MATLAB_MODULE=matlab

# How run_all.sh chains the jobs (preproc -> prep -> GLM):
#   afterok   next stage runs only if ALL tasks succeeded (safest)
#   afterany  next stage runs regardless of exit status
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
# subject-sessions at runtime. Rule of thumb: GRIDCAT_CPUS minus 1.
MAX_WORKERS=29

# Stage data to node-local scratch to bypass NFS? Input is copied to the
# node's /tmp, output written locally, results copied back at the end.
USE_LOCAL_SCRATCH=true

# Local tmp to request from SLURM (--gres=tmp:SIZE). Must exceed your
# GLM_runauto directory (~13G) plus expected output. Empty skips the request.
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
