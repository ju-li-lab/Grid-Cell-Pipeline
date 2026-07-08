#!/usr/bin/env bash
# ==============================================================================
# ASHS Label Extraction — Configuration
# ==============================================================================
# Edit this file to control what gets extracted and where it goes.
# Then run:  bash ashs_extract_labels.sh
#
# You do NOT need coding experience to use this. Just fill in the values below.
# Lines starting with # are comments and are ignored.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1) INPUT: Where is the ASHS output?
# ------------------------------------------------------------------------------
# Full path to the top-level ASHS output directory.
# This directory should contain sub-XX folders.
#
# Example:  /sc-projects/sc-proj-cc02-brainspace/sirius_wip/ashs_output
ASHS_DIR="/charite-store-f/f-cc04-stahn_lab/projects/BRAVITY/Data/BRACE/mri/bids/derivatives/ashs/abc_prisma"

# ------------------------------------------------------------------------------
# 2) WHICH SEGMENTATION TYPE to use?
# ------------------------------------------------------------------------------
# ASHS produces three segmentation variants per hemisphere. Pick ONE:
#
#   corr_usegray   — corrected with gray matter info (RECOMMENDED for most uses)
#   corr_nogray    — corrected without gray matter info
#   heur           — heuristic-based segmentation
#
SEG_TYPE="corr_usegray"

# ------------------------------------------------------------------------------
# 3) WHICH LABELS to extract?
# ------------------------------------------------------------------------------
# List the label numbers you want extracted, separated by spaces.
# Each label becomes its own binary mask file (1 = region, 0 = everything else).
#
# Label reference for this atlas:
#   0  = Background (do not extract)
#   1  = CA1
#   2  = CA2
#   3  = DG (Dentate Gyrus)
#   4  = CA3
#   5  = Tail (Hippocampal Tail)
#   6  = Hippocampal_sulcus
#   7  = MISC
#   8  = SUB (Subiculum)
#   9  = ERC (Entorhinal Cortex)
#   10 = BA35
#   11 = BA36
#   12 = PHC (Parahippocampal Cortex)
#   13 = MISC2
#   14 = CS (Collateral Sulcus)
#
# Example: extract ERC only
LABELS_TO_EXTRACT="9"
#
# Example: extract ERC + PHC + SUB
# LABELS_TO_EXTRACT="9 12 8"
#
# Example: extract all hippocampal subfields
# LABELS_TO_EXTRACT="1 2 3 4 5"

# Friendly names for the labels (must match order of LABELS_TO_EXTRACT).
# These are used in output filenames. Use simple names, no spaces.
# If left empty, the script uses the default names from the atlas.
LABEL_NAMES="ErC"
# Example: LABEL_NAMES="ERC PHC SUB"

# ------------------------------------------------------------------------------
# 4) WHICH HEMISPHERES?
# ------------------------------------------------------------------------------
# Which hemispheres to process:  left  right  both
HEMISPHERES="both"

# ------------------------------------------------------------------------------
# 5) OUTPUT MODE: flat folder or BIDS directory?
# ------------------------------------------------------------------------------
# Choose ONE output mode:
#
#   flat   — All extracted masks go into a single flat output folder.
#            Good for quick access and GridCAT input.
#
#   bids   — Extracted masks are copied into the correct
#            sub-XX/ses-YY/anat/ folders of a BIDS directory.
#            Good for keeping things organized with your main dataset.
#
OUTPUT_MODE="bids"

# For OUTPUT_MODE="flat":
# Where to put all extracted masks. Will be created if it does not exist.
FLAT_OUTPUT_DIR="/sc-projects/sc-proj-cc02-brainspace/myspace/julius/BRACE_GC_testing/BRAVE/extracted_masks"

# For OUTPUT_MODE="bids":
# Path to the BIDS dataset root (must already contain sub-XX/ses-YY/anat/).
# Masks will be placed into the matching sub/ses/anat folder.
BIDS_OUTPUT_DIR="/sc-projects/sc-proj-cc02-brainspace/myspace/julius/BRACE_GC_testing/BRACE/bids"

# ------------------------------------------------------------------------------
# 6) MULTIPLE RUNS — how to handle them
# ------------------------------------------------------------------------------
# When ASHS was run with multiple T1/T2 pairings for the same subject+session,
# the output contains multiple run directories (e.g., t1-run01_t2-run01).
#
# Options:
#   ask      — The script will pause and ask you which run to use (RECOMMENDED)
#   first    — Automatically pick the first run found (alphabetical order)
#   all      — Extract from ALL runs (output filenames will include run info)
#
MULTI_RUN_MODE="all"

# ------------------------------------------------------------------------------
# 7) SUBJECT/SESSION FILTER (optional)
# ------------------------------------------------------------------------------
# Leave empty to process ALL subjects and sessions found in ASHS_DIR.
# Or list specific subjects/sessions to process.
#
# Examples:
#   SUBJECTS="sub-01c06 sub-02c06"       # only these subjects
#   SESSIONS="ses-01"                      # only session 01
#   SUBJECTS=""                            # all subjects (default)
#   SESSIONS=""                            # all sessions (default)
#
SUBJECTS=""
SESSIONS=""

# ------------------------------------------------------------------------------
# 8) ADVANCED OPTIONS
# ------------------------------------------------------------------------------
# Combine extracted labels into one multi-label mask?
# If "yes", in addition to individual binary masks, a combined mask file
# is created that contains all requested labels in one file.
COMBINE_LABELS="no"

# Output file format:  nii.gz  or  nii
OUTPUT_FORMAT="nii"

# Overwrite existing output files?  yes / no
OVERWRITE="no"

# Enable dry-run mode? (shows what would happen without actually doing it)
DRY_RUN="no"
