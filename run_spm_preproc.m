function run_spm_preproc(BIDS_ROOT, SUB, SES, cfgFile, stages)
% run_spm_preproc(BIDS_ROOT, SUB, SES, cfgFile, stages)
%
% Example:
%   run_spm_preproc('/sc-projects/.../b2_bids','sub-01s06','ses-01')
%   run_spm_preproc('/sc-projects/.../b2_bids','sub-01s06','ses-01', '/path/to/pipeline_config.cfg')
%   run_spm_preproc('/sc-projects/.../b2_bids','sub-01s06','ses-01', '', 'realign,coreg,smooth')
%
% For one subject/session:
% - Picks one run of each input (BOLD per task, T2w, fieldmap, reverse-PE EPI)
%     and copies it into the derivatives directory
% - Finds fmap magnitude1 + phasediff, or the reverse phase-encode EPI
% - For each task, builds a voxel displacement map matched to that run's geometry
% - Runs Realign & Unwarp for the whole session with one data block per task
%     (or Realign only, depending on PREPROC_MODE in config)
% - Creates a session mean by averaging the per-task means
% - Coregisters T2w -> session mean and reslices ROI left/right into EPI space
% - Optionally smooths the preprocessed BOLD files
%
% WHERE THE OUTPUT GOES
% ---------------------
% Nothing is ever written into BIDS_ROOT. Every input SPM touches is copied
% into a BIDS-style derivatives tree first:
%
%   <DERIV_ROOT>/spm-preproc/<sub>/<ses>/func/   u*_bold.nii, rp_*.txt, vdm_*,
%                                                topup_*, meanu_session.nii, su*
%                                     /anat/     the T2w and the resliced ROIs
%                                     /fmap/     the fieldmap inputs
%
% This matters beyond tidiness: SPM's Realign & Unwarp and Coregister write the
% transforms they estimate into the *headers of their input images*, and
% ensure_nii used to leave decompressed twins of every .nii.gz behind. Reading
% raw BIDS and writing derivatives keeps the source data untouched.
%
% A full run (a stage list that starts the preprocessing from its beginning)
% clears the session's derivatives first, so the result is never a mixture of
% this run and the last one. Partial runs keep them — that is what they resume
% from. See DERIV_RESET in pipeline_config.cfg.
%
% MULTIPLE RUNS OF THE SAME SCAN
% ------------------------------
% When a session holds several runs of a scan, one is chosen per input, in this
% order of preference:
%
%   1. The choice in run_selection.tsv (RUN_SELECTION_FILE). Write it with
%      bash scan_multirun.sh, which lists every run with its acquisition time.
%   2. The run acquired closest in time to the BOLD being preprocessed, when
%      the JSON sidecars carry AcquisitionTime. Run numbers are counted per
%      modality, so run-2 of the T2w need not belong with run-2 of a task —
%      if the subject climbed out of the scanner in between, they do not.
%      Acquisition time is what actually says which scans belong together.
%   3. Nothing else is safe. With STRICT_RUN_MATCHING=true the run stops and
%      asks for a selection; otherwise it takes the highest run number and
%      says loudly that it guessed.
%
% The chosen runs are recorded in <ses>/preproc_run_selection.tsv, which the
% GridCAT preparation reads so the event files line up with the right run.
%
% STAGES (5th argument, or PREPROC_STAGES in pipeline_config.cfg)
% ---------------------------------------------------------------
% The work above is split into five stages that can be run independently. Every
% stage reads what it needs from disk, so you can stop after any stage, inspect
% or replace the intermediate files, and resume later.
%
%   topup    Estimate the distortion field with FSL topup. Writes
%            <func>/topup_results_* and <func>/topup_acqparams.txt.
%            Does nothing unless PREPROC_MODE=topup.
%   vdm      Build one voxel displacement map per task: <func>/vdm_task-<N>.nii.
%            Built from the GRE fieldmap (PREPROC_MODE=realign_unwarp) or from
%            the topup field (PREPROC_MODE=topup, TOPUP_APPLY_METHOD=vdm).
%            Does nothing for realign_only or for topup+applytopup.
%   realign  Apply the correction and do motion correction: Realign & Unwarp,
%            or applytopup followed by Realign, or Realign alone. Writes
%            u*_bold.nii, rp_*.txt, the per-task means and meanu_session.nii.
%   coreg    Coregister T2w -> session mean and reslice the ROI masks into EPI
%            space (r*_mask.nii). Does nothing when DO_COREG_ROIS=0.
%   smooth   Gaussian smoothing of the preprocessed BOLD files (su*_bold.nii).
%            Does nothing when SMOOTH_FWHM=0.
%
% Accepted values: a comma-separated list of stage names, or one of the
% shorthands 'all' (default), 'fieldmap' (= topup,vdm) and 'post_fieldmap'
% (= realign,coreg,smooth). Order does not matter — stages always run in the
% canonical order listed above.
%
% RESUMING FROM FIELDMAPS YOU MADE YOURSELF
% -----------------------------------------
% Put your FSL topup output on the cluster and point TOPUP_EXISTING_PREFIX at it
% (or copy it into <func>/ as topup_results_*), then run the stages
% 'vdm,realign,coreg,smooth'. Alternatively drop ready-made voxel displacement
% maps in as <func>/vdm_task-<N>.nii and run 'realign,coreg,smooth'.
%
% All settings come from pipeline_config.cfg via read_pipeline_config.m

% --------- Handle optional arguments ----------
if nargin < 4, cfgFile = ''; end
if nargin < 5, stages  = ''; end

if isempty(cfgFile)
    % Auto-detect pipeline_config.cfg in same directory as this script
    scriptDir = fileparts(mfilename('fullpath'));
    cfgFile = fullfile(scriptDir, 'pipeline_config.cfg');
end

assert(isfile(cfgFile), 'pipeline_config.cfg not found at: %s\nProvide cfgFile as 4th argument to run_spm_preproc.', cfgFile);

% --------- Read configuration ----------
cfgRaw = read_pipeline_config(cfgFile);

% --------- Map config keys to CFG struct ----------
CFG = struct();

% Paths
CFG.SPM_DIR = cfgRaw.SPM_DIR;

% Acquisition parameters (note: config stores in MS, some script vars may differ)
CFG.TR = cfgRaw.TR;                            % seconds
CFG.TOTAL_READOUT = cfgRaw.TOTAL_READOUT_MS;  % milliseconds
CFG.BLIPDIR = cfgRaw.BLIP_DIRECTION;
CFG.EPI_BASED_FIELDMAP = cfgRaw.EPI_BASED_FIELDMAP;
CFG.TE_SHORT = cfgRaw.TE_SHORT_MS;             % milliseconds
CFG.TE_LONG = cfgRaw.TE_LONG_MS;               % milliseconds
CFG.MASKBRAIN = cfgRaw.VDM_MASKBRAIN;

% Tasks — can be numeric array [1 2 3] or string 'run1,run2'
% Normalise to cell array of strings: {'run1','run2'} or {'1','2','3'}
if isnumeric(cfgRaw.TASKS)
    CFG.TASK_LABELS = arrayfun(@num2str, cfgRaw.TASKS, 'UniformOutput', false);
else
    CFG.TASK_LABELS = strtrim(strsplit(cfgRaw.TASKS, ','));
end

% Run selection file (for multi-run T2w/BOLD/fieldmap/reverse-PE)
if isfield(cfgRaw, 'RUN_SELECTION_FILE') && ~isempty(cfgRaw.RUN_SELECTION_FILE)
    CFG.RUN_SELECTION_FILE = cfgRaw.RUN_SELECTION_FILE;
else
    CFG.RUN_SELECTION_FILE = '';
end

% How far apart two scans may be acquired and still count as the same visit to
% the scanner. Runs further apart than this are treated as different blocks,
% which is what happens when the subject gets out and comes back.
if isfield(cfgRaw, 'RUN_MATCH_GAP_MIN') && ~isempty(cfgRaw.RUN_MATCH_GAP_MIN)
    CFG.RUN_MATCH_GAP_SEC = double(cfgRaw.RUN_MATCH_GAP_MIN) * 60;
else
    CFG.RUN_MATCH_GAP_SEC = 20 * 60;
end

% Stop rather than guess when several runs exist and nothing decides between
% them. Silently pairing the wrong T2w with a task is worse than not running.
CFG.STRICT_RUN_MATCHING = cfg_flag(cfgRaw, 'STRICT_RUN_MATCHING', true);

% Clearing the derivatives before a full run: auto | always | never
if isfield(cfgRaw, 'DERIV_RESET') && ~isempty(cfgRaw.DERIV_RESET)
    CFG.DERIV_RESET = lower(strtrim(char(cfgRaw.DERIV_RESET)));
else
    CFG.DERIV_RESET = 'auto';
end

% ROI patterns
CFG.ROI_PATTERN_LEFT = cfgRaw.ROI_PATTERN_LEFT;
CFG.ROI_PATTERN_RIGHT = cfgRaw.ROI_PATTERN_RIGHT;
CFG.DO_COREG_ROIS = cfgRaw.DO_COREG_ROIS;

% Preprocessing mode
CFG.PREPROC_MODE = cfgRaw.PREPROC_MODE;  % 'realign_unwarp', 'topup', or 'realign_only'

% FSL topup parameters (only used when PREPROC_MODE=topup)
if strcmp(CFG.PREPROC_MODE, 'topup')
    CFG.TOPUP_PE_DIR_BOLD = cfgRaw.TOPUP_PE_DIR_BOLD;
    CFG.TOPUP_PE_DIR_REVERSE = cfgRaw.TOPUP_PE_DIR_REVERSE;
    CFG.TOPUP_READOUT_SEC = cfgRaw.TOPUP_READOUT_SEC;
    CFG.TOPUP_REVERSE_PE_PATTERN = cfgRaw.TOPUP_REVERSE_PE_PATTERN;
    if isfield(cfgRaw, 'TOPUP_REVERSE_PE_DIR') && ~isempty(cfgRaw.TOPUP_REVERSE_PE_DIR)
        CFG.TOPUP_REVERSE_PE_DIR = cfgRaw.TOPUP_REVERSE_PE_DIR;
    else
        CFG.TOPUP_REVERSE_PE_DIR = 'auto';
    end
    if isfield(cfgRaw, 'TOPUP_CONFIG') && ~isempty(cfgRaw.TOPUP_CONFIG)
        CFG.TOPUP_CONFIG = cfgRaw.TOPUP_CONFIG;
    else
        CFG.TOPUP_CONFIG = '';
    end
    CFG.TOPUP_APPLY_METHOD = cfgRaw.TOPUP_APPLY_METHOD;
    if isfield(cfgRaw, 'TOPUP_INTERP') && ~isempty(cfgRaw.TOPUP_INTERP)
        CFG.TOPUP_INTERP = cfgRaw.TOPUP_INTERP;
    else
        CFG.TOPUP_INTERP = 'spline';
    end
end

% Precalculated fieldmap parameters (only used when PREPROC_MODE=precalc_fieldmap)
% The fieldmap and its magnitude are what make_fieldmaps.sh writes into fmap/.
% The extension is appended by the search, so strip it if the pattern was
% written with one (the ROI patterns in the same config do include it).
stripExt = @(x) regexprep(x, '\.nii(\.gz)?$', '');
if isfield(cfgRaw, 'FIELDMAP_PATTERN') && ~isempty(cfgRaw.FIELDMAP_PATTERN)
    CFG.FIELDMAP_PATTERN = stripExt(cfgRaw.FIELDMAP_PATTERN);
else
    CFG.FIELDMAP_PATTERN = '_fieldmap';
end
if isfield(cfgRaw, 'FIELDMAP_MAGNITUDE_PATTERN') && ~isempty(cfgRaw.FIELDMAP_MAGNITUDE_PATTERN)
    CFG.FIELDMAP_MAGNITUDE_PATTERN = stripExt(cfgRaw.FIELDMAP_MAGNITUDE_PATTERN);
else
    CFG.FIELDMAP_MAGNITUDE_PATTERN = '_magnitude';
end
% Units of the fieldmap on disk. fsl_prepare_fieldmap writes rad/s; SPM's
% FieldMap toolbox expects Hz, so rad/s is divided by 2*pi before use.
if isfield(cfgRaw, 'FIELDMAP_UNITS') && ~isempty(cfgRaw.FIELDMAP_UNITS)
    CFG.FIELDMAP_UNITS = lower(strtrim(cfgRaw.FIELDMAP_UNITS));
else
    CFG.FIELDMAP_UNITS = 'rad/s';
end

% Session mean settings
CFG.DO_SESSION_MEAN = 1;
CFG.SESSION_MEAN_NAME = 'meanu_session.nii';

% VDM / FieldMap parameters
CFG.VDM_WRAP = cfgRaw.VDM_WRAP;           % Should be [0 1 0] from config
CFG.VDM_MASK = cfgRaw.VDM_MASK;
CFG.VDMFLAGS_BASE = struct('udir', 1, 'rinterp', 4, 'wrap', CFG.VDM_WRAP, ...
                           'mask', CFG.VDM_MASK, 'order', 1, 'prefix', 'vdm_');

% Reslice / Unwarp options
CFG.WHICHUNWARP = cfgRaw.RESLICE_WHICH;   % Should be [2 1] from config
CFG.RESLICE_INTERP = cfgRaw.RESLICE_INTERP;
CFG.RESLICE_WRAP = cfgRaw.RESLICE_WRAP;
CFG.RESLICE_MASK = cfgRaw.RESLICE_MASK;
CFG.RESLICE_PREFIX = cfgRaw.RESLICE_PREFIX;

% Realign & Unwarp estimation options
CFG.REALIGN_QUALITY = cfgRaw.REALIGN_QUALITY;
CFG.REALIGN_SEP = cfgRaw.REALIGN_SEP;
CFG.REALIGN_FWHM = cfgRaw.REALIGN_FWHM;
CFG.REALIGN_RTM = cfgRaw.REALIGN_RTM;
CFG.REALIGN_EINTERP = cfgRaw.REALIGN_EINTERP;
CFG.REALIGN_EWRAP = cfgRaw.REALIGN_EWRAP;

% Unwarp estimation options
CFG.UNWARP_BASFCN = cfgRaw.UNWARP_BASFCN;
CFG.UNWARP_REGORDER = cfgRaw.UNWARP_REGORDER;
CFG.UNWARP_LAMBDA = cfgRaw.UNWARP_LAMBDA;
CFG.UNWARP_JM = cfgRaw.UNWARP_JM;
CFG.UNWARP_FOT = cfgRaw.UNWARP_FOT;
CFG.UNWARP_SOT = cfgRaw.UNWARP_SOT;
CFG.UNWARP_UWFWHM = cfgRaw.UNWARP_UWFWHM;
CFG.UNWARP_REM = cfgRaw.UNWARP_REM;
CFG.UNWARP_NOI = cfgRaw.UNWARP_NOI;
CFG.UNWARP_EXPROUND = cfgRaw.UNWARP_EXPROUND;

% VDM unwarp flags
CFG.VDM_UFLAGS_METHOD = cfgRaw.VDM_UFLAGS_METHOD;
CFG.VDM_UFLAGS_FWHM = cfgRaw.VDM_UFLAGS_FWHM;
CFG.VDM_UFLAGS_PAD = cfgRaw.VDM_UFLAGS_PAD;
CFG.VDM_UFLAGS_WS = cfgRaw.VDM_UFLAGS_WS;

% Coregistration options
CFG.COREG_COST_FUN = cfgRaw.COREG_COST_FUN;
CFG.COREG_SEP = cfgRaw.COREG_SEP;
CFG.COREG_FWHM = cfgRaw.COREG_FWHM;
CFG.COREG_ROI_INTERP = cfgRaw.COREG_ROI_INTERP;

% Smoothing (optional — off by default)
if isfield(cfgRaw, 'SMOOTH_FWHM') && ~isempty(cfgRaw.SMOOTH_FWHM)
    CFG.SMOOTH_FWHM = cfgRaw.SMOOTH_FWHM;
else
    CFG.SMOOTH_FWHM = 0;  % default: no smoothing
end
if isfield(cfgRaw, 'SMOOTH_PREFIX') && ~isempty(cfgRaw.SMOOTH_PREFIX)
    CFG.SMOOTH_PREFIX = cfgRaw.SMOOTH_PREFIX;
else
    CFG.SMOOTH_PREFIX = 's';
end

% Pre-computed FSL topup output (optional).
% Points at topup results you produced yourself, so the "topup" stage can be
% skipped. {SUB} and {SES} are substituted. Empty = use <func>/topup_results.
if isfield(cfgRaw, 'TOPUP_EXISTING_PREFIX') && ~isempty(cfgRaw.TOPUP_EXISTING_PREFIX)
    CFG.TOPUP_EXISTING_PREFIX = cfgRaw.TOPUP_EXISTING_PREFIX;
else
    CFG.TOPUP_EXISTING_PREFIX = '';
end

% --------- Resolve which stages to run ----------
% Explicit argument wins over PREPROC_STAGES in the config file.
if isempty(stages) && isfield(cfgRaw, 'PREPROC_STAGES')
    stages = cfgRaw.PREPROC_STAGES;
end
STAGES = resolve_stages(stages);

% --------- Safety / init ----------
assert(isfolder(BIDS_ROOT), 'BIDS_ROOT not found: %s', BIDS_ROOT);
assert(startsWith(SUB,'sub-') && startsWith(SES,'ses-'), 'SUB/SES must look like sub-01s06 / ses-01');
assert(ismember(CFG.PREPROC_MODE, {'realign_unwarp','topup','realign_only','precalc_fieldmap'}), ...
    'Unknown PREPROC_MODE: %s (must be realign_unwarp, topup, precalc_fieldmap or realign_only)', ...
    CFG.PREPROC_MODE);
if strcmp(CFG.PREPROC_MODE, 'precalc_fieldmap')
    assert(ismember(CFG.FIELDMAP_UNITS, {'rad/s','rads','rad_per_s','hz'}), ...
        'Unknown FIELDMAP_UNITS: %s (must be rad/s or Hz)', CFG.FIELDMAP_UNITS);
end
if strcmp(CFG.PREPROC_MODE, 'topup')
    assert(ismember(CFG.TOPUP_APPLY_METHOD, {'applytopup','vdm'}), ...
        'Unknown TOPUP_APPLY_METHOD: %s (must be applytopup or vdm)', CFG.TOPUP_APPLY_METHOD);
end

addpath(CFG.SPM_DIR);
spm('defaults','FMRI');
spm_jobman('initcfg');

% --------- Raw BIDS: read only, never written to ----------
RAW = struct();
RAW.ses  = fullfile(BIDS_ROOT, SUB, SES);
RAW.func = fullfile(RAW.ses, 'func');
RAW.anat = fullfile(RAW.ses, 'anat');
RAW.fmap = fullfile(RAW.ses, 'fmap');

assert(isfolder(RAW.func), 'Missing func dir: %s', RAW.func);

% --------- Derivatives: everything this run writes ----------
DER = bids_deriv(cfgRaw, SUB, SES);
doReset = deriv_reset_wanted(CFG, STAGES);
bids_reset_deriv(DER, doReset);

funcDir = DER.func;   % kept as short names: every stage below writes here
anatDir = DER.anat;
fmapDir = DER.fmap;

fprintf('\n=== %s / %s ===\n', SUB, SES);
fprintf('Preprocessing mode: %s\n', CFG.PREPROC_MODE);
fprintf('Stages to run:      %s\n', strjoin(STAGES.list, ' -> '));
if ~isempty(STAGES.skipped)
    fprintf('Stages skipped:     %s\n', strjoin(STAGES.skipped, ', '));
end
fprintf('Reading from:       %s\n', RAW.ses);
fprintf('Writing to:         %s%s\n', DER.sesDir, tern(doReset, '   (cleared first)', '   (resuming)'));

% --------- Which directories a run actually needs depends on the stages ------
if STAGES.vdm && ismember(CFG.PREPROC_MODE, {'realign_unwarp','precalc_fieldmap'})
    assert(isfolder(RAW.fmap) || isfolder(DER.fieldmapSes), ...
        ['No fieldmap directory for %s / %s.\n' ...
         '  Looked in %s\n  and in %s'], SUB, SES, RAW.fmap, DER.fieldmapSes);
elseif STAGES.topup && strcmp(CFG.PREPROC_MODE, 'topup') && strcmp(CFG.TOPUP_REVERSE_PE_DIR, 'fmap')
    assert(isfolder(RAW.fmap), 'Missing fmap dir: %s (TOPUP_REVERSE_PE_DIR=fmap)', RAW.fmap);
end
if STAGES.coreg && CFG.DO_COREG_ROIS
    assert(isfolder(RAW.anat), 'Missing anat dir: %s', RAW.anat);
end

% --------- Load run selection (multi-run overrides) ----------
runSel = bids_load_selection(CFG.RUN_SELECTION_FILE, SUB, SES);
if ~isempty(CFG.RUN_SELECTION_FILE) && isfile(CFG.RUN_SELECTION_FILE)
    fprintf('Run selection:      %s\n', CFG.RUN_SELECTION_FILE);
end

% =====================================================================
%  Choose one run of each input, and copy it into the derivatives
% =====================================================================
nTasks = numel(CFG.TASK_LABELS);
taskVols      = cell(nTasks,1);
taskVDM       = cell(nTasks,1);
taskBoldFiles = cell(nTasks,1);

needVols = STAGES.vdm || STAGES.realign;

% Smoothing on its own reads only what the earlier stages left in the
% derivatives, so there is nothing to choose and nothing to copy.
needSelection = STAGES.topup || STAGES.vdm || STAGES.realign || STAGES.coreg;

boldPicks = repmat(struct('label', '', 'path', '', 'runLabel', '', ...
                          'acqSec', NaN, 'reason', ''), nTasks, 1);
anchor = struct('sec', NaN, 'secs', NaN, 'name', 'the functional data');
for ii = 1:nTasks
    boldPicks(ii).label = CFG.TASK_LABELS{ii};
    taskVDM{ii} = '';
end

if needSelection
    % --------- The BOLD run per task: what everything else matches to -------
    fprintf('\n--- Choosing runs ---\n');
    [boldPicks, anchor] = select_task_bolds(RAW.func, SUB, SES, CFG, runSel, needVols);

    for ii = 1:nTasks
        if isempty(boldPicks(ii).path), continue; end
        taskBoldFiles{ii} = bids_stage(boldPicks(ii).path, funcDir);
    end

    if ~isnan(anchor.sec)
        fprintf('  Matching other scans to %s (acquired %s)\n', anchor.name, clock_of(anchor.sec));
    else
        fprintf('  No acquisition times in the BOLD sidecars — other scans cannot be matched by time.\n');
    end
end

% --------- Fieldmap inputs (only for the stages that consume them) -----------
phasemap = '';
magnitude1 = '';
reversePE_epi = '';
fieldmapFile = '';
fieldmapMag = '';
fieldmapSrc = '';       % raw path of the chosen fieldmap, for the manifest
reverseSrc  = '';
phasediffSrc = '';

if strcmp(CFG.PREPROC_MODE, 'realign_unwarp') && STAGES.vdm
    % GRE fieldmap mode: need phasediff + the magnitude it pairs with
    [phasediffSrc, magSrc] = select_phasediff(RAW.fmap, SUB, SES, CFG, runSel, anchor);
    phasemap   = bids_stage(phasediffSrc, fmapDir);
    magnitude1 = bids_stage(magSrc,       fmapDir);

    fprintf('Fieldmap phasediff:  %s\n', phasemap);
    fprintf('Fieldmap magnitude1: %s\n', magnitude1);

elseif strcmp(CFG.PREPROC_MODE, 'precalc_fieldmap') && STAGES.vdm
    % Precalculated fieldmap mode: a ready-made B0 map plus the magnitude image
    % it was unwrapped against — what make_fieldmaps.sh writes.
    [fieldmapSrc, magSrc] = select_precalc_fieldmap(RAW.fmap, DER.fieldmapSes, SUB, SES, CFG, runSel, anchor);
    fieldmapFile = bids_stage(fieldmapSrc, fmapDir);
    fieldmapMag  = bids_stage(magSrc,      fmapDir);
    fprintf('Fieldmap:           %s\n', fieldmapFile);
    fprintf('Fieldmap magnitude: %s\n', fieldmapMag);
    fprintf('Fieldmap units:     %s\n', CFG.FIELDMAP_UNITS);

elseif strcmp(CFG.PREPROC_MODE, 'topup') && STAGES.topup
    % Topup mode: find the reverse-PE EPI acquired with these functional runs.
    reverseSrc = select_reverse_pe(RAW, SUB, SES, CFG, runSel, anchor);
    reversePE_epi = bids_stage(reverseSrc, fmapDir);
    fprintf('Reverse-PE EPI:  %s\n', reversePE_epi);
    fprintf('Topup apply method: %s\n', CFG.TOPUP_APPLY_METHOD);
end

% --------- T2w + ROIs (only the coreg stage uses them) ----------
t2w = '';  t2wSrc = '';
roiL = ''; roiR = ''; hasRoiL = false; hasRoiR = false;
roiLsrc = ''; roiRsrc = '';

if STAGES.coreg && CFG.DO_COREG_ROIS
    t2wSrc = select_t2w(RAW.anat, SUB, SES, CFG, runSel, anchor);
    t2w = bids_stage(t2wSrc, anatDir);
    fprintf('T2w: %s\n', t2w);

    % ROI masks may have been imported into the derivatives by move_rois.sh, or
    % may still be sitting in the raw anat/ from an earlier way of working.
    [roiLsrc, roiRsrc, hasRoiL, hasRoiR] = find_rois({DER.roiSes, RAW.anat}, SUB, SES, CFG, runSel, anchor);
    if hasRoiL, roiL = bids_stage(roiLsrc, anatDir); end
    if hasRoiR, roiR = bids_stage(roiRsrc, anatDir); end
    fprintf('ROI left exists:  %d', hasRoiL);
    if hasRoiL, fprintf(' (%s)', roiL); end
    fprintf('\n');
    fprintf('ROI right exists: %d', hasRoiR);
    if hasRoiR, fprintf(' (%s)', roiR); end
    fprintf('\n');
end

% --------- Record what was chosen ----------
% The GridCAT preparation reads this back, so the event tables end up with the
% run that was actually preprocessed even when their filenames say nothing
% about runs. A partial run only decided some of these, so the file is merged
% rather than replaced — otherwise re-running the coreg stage on its own would
% erase which BOLD run the realignment used.
if needSelection
    write_run_manifest(DER, SUB, SES, CFG, boldPicks, ...
        struct('T2w', t2wSrc, 'fieldmap', fieldmapSrc, 'phasediff', phasediffSrc, ...
               'reverse', reverseSrc, 'roi_left', roiLsrc, 'roi_right', roiRsrc), ...
        anchor);
end

% --------- Validate config against the BIDS JSON sidecars we chose ----------
% Reading the sidecars of the selected files (not "any file in the folder")
% means the check describes the run actually being preprocessed.
if STAGES.topup || STAGES.vdm || STAGES.realign
    validate_bids_params(taskBoldFiles, reversePE_epi, phasemap, fieldmapFile, CFG);
end

% --------- Expand the staged BOLD into per-volume references ----------
for ii = 1:nTasks
    tLabel = CFG.TASK_LABELS{ii};
    fprintf('\n--- Task %s ---\n', tLabel);
    if isempty(taskBoldFiles{ii})
        fprintf('BOLD: (not needed by the selected stages)\n');
        continue;
    end
    fprintf('BOLD: %s   [%s]\n', taskBoldFiles{ii}, boldPicks(ii).reason);

    if needVols
        taskVols{ii} = expand_4d(taskBoldFiles{ii});   % '...nii,1' '...nii,2' ...
    end
end

% Does this configuration need voxel displacement maps at all?
needVDM = ismember(CFG.PREPROC_MODE, {'realign_unwarp','precalc_fieldmap'}) || ...
          (strcmp(CFG.PREPROC_MODE, 'topup') && strcmp(CFG.TOPUP_APPLY_METHOD, 'vdm'));

% =====================================================================
%  STAGE topup — estimate the distortion field with FSL
% =====================================================================
topupPrefix  = '';
topupAcqFile = '';

if strcmp(CFG.PREPROC_MODE, 'topup') && (STAGES.topup || STAGES.vdm || STAGES.realign)
    topupPrefix = resolve_topup_prefix(funcDir, SUB, SES, CFG, ~STAGES.topup);

    if STAGES.topup
        fprintf('\n=== STAGE topup: FSL topup field estimation ===\n');
        run_fsl_topup(taskBoldFiles{1}, reversePE_epi, topupPrefix, CFG);
    else
        fprintf('\n=== STAGE topup: skipped, reusing existing topup output ===\n');
        assert_topup_outputs(topupPrefix, CFG, STAGES);
    end
    fprintf('Topup prefix: %s\n', topupPrefix);

    if STAGES.realign && strcmp(CFG.TOPUP_APPLY_METHOD, 'applytopup')
        topupAcqFile = resolve_acqparams(topupPrefix, funcDir, CFG);
        fprintf('Topup acqparams: %s\n', topupAcqFile);
    end
end

% =====================================================================
%  STAGE vdm — one voxel displacement map per task
% =====================================================================
if STAGES.vdm
    if ~needVDM
        fprintf('\n=== STAGE vdm: nothing to build in mode %s ===\n', CFG.PREPROC_MODE);
    else
        fprintf('\n=== STAGE vdm: building voxel displacement maps ===\n');

        % SPM's FieldMap toolbox works in Hz. fsl_prepare_fieldmap writes rad/s,
        % so convert once per session before looping over the tasks.
        fieldmapHz = '';
        if strcmp(CFG.PREPROC_MODE, 'precalc_fieldmap')
            fieldmapHz = fieldmap_to_hz(fieldmapFile, funcDir, CFG);
        end

        for ii = 1:nTasks
            tLabel = CFG.TASK_LABELS{ii};
            epiRef = taskVols{ii}{1};   % first volume as VDM reference
            switch CFG.PREPROC_MODE
                case 'realign_unwarp'
                    fprintf('Calculating VDM for task-%s...\n', tLabel);
                    taskVDM{ii} = calc_vdm_for_run_to_func(phasemap, magnitude1, epiRef, funcDir, fmapDir, CFG, tLabel);
                case 'precalc_fieldmap'
                    fprintf('Calculating VDM from precalculated fieldmap for task-%s...\n', tLabel);
                    taskVDM{ii} = calc_vdm_from_fieldmap(fieldmapHz, fieldmapMag, epiRef, funcDir, fmapDir, CFG, tLabel);
                otherwise
                    fprintf('Converting topup field to VDM for task-%s...\n', tLabel);
                    taskVDM{ii} = convert_topup_to_vdm(topupPrefix, epiRef, funcDir, CFG, tLabel);
            end
            fprintf('VDM (task-%s): %s\n', tLabel, taskVDM{ii});
        end

        % Drop the Hz intermediate; the VDMs are what the realign stage reads
        if ~isempty(fieldmapHz) && isfile(fieldmapHz)
            delete(fieldmapHz);
        end
    end
end

% =====================================================================
%  STAGE realign — apply the correction + motion correction
% =====================================================================
sessionMean     = '';
haveSessionMean = false;

if STAGES.realign
    fprintf('\n=== STAGE realign: %s ===\n', CFG.PREPROC_MODE);

    if needVDM
        % Fill in any VDM this run did not build itself (vdm stage skipped)
        for ii = 1:nTasks
            if isempty(taskVDM{ii})
                taskVDM{ii} = existing_task_vdm(funcDir, CFG.TASK_LABELS{ii}, ii);
                assert(~isempty(taskVDM{ii}), [ ...
                    'No voxel displacement map for task-%s.\n' ...
                    'Expected: %s\n' ...
                    'Run the "vdm" stage first, or copy a ready-made VDM to that path.'], ...
                    CFG.TASK_LABELS{ii}, fullfile(funcDir, sprintf('vdm_task-%s.nii', CFG.TASK_LABELS{ii})));
                fprintf('Reusing existing VDM (task-%s): %s\n', CFG.TASK_LABELS{ii}, taskVDM{ii});
            end
        end
    end

    if strcmp(CFG.PREPROC_MODE, 'topup') && strcmp(CFG.TOPUP_APPLY_METHOD, 'applytopup')
        % Correct each task's 4D BOLD with applytopup, then realign only
        fprintf('\n--- Applying topup correction via applytopup ---\n');
        correctedVols = cell(nTasks,1);
        for ii = 1:nTasks
            tLabel = CFG.TASK_LABELS{ii};
            % Keep every BIDS entity of the file it came from, so a multi-run
            % session still says which run was corrected.
            correctedBold = regexprep(taskBoldFiles{ii}, '_bold\.nii$', '_bold_dc.nii');
            correctedBold = run_fsl_applytopup(taskBoldFiles{ii}, topupPrefix, topupAcqFile, correctedBold, 1, CFG);
            fprintf('Corrected BOLD (task-%s): %s\n', tLabel, correctedBold);
            correctedVols{ii} = expand_4d(correctedBold);
        end

        % Realign only (distortion already corrected by applytopup)
        run_realign_only_multi(correctedVols, CFG);

    elseif needVDM
        run_realign_unwarp_multi(taskVols, taskVDM, CFG);

    else
        % realign_only mode
        run_realign_only_multi(taskVols, CFG);
    end

    % Per-task means -> session mean (the coregistration reference)
    sessionMean = resolve_session_mean(funcDir, CFG, true);
    haveSessionMean = true;
end

% =====================================================================
%  STAGE coreg — T2w -> session mean, reslice ROI masks into EPI space
% =====================================================================
if STAGES.coreg
    if ~CFG.DO_COREG_ROIS
        fprintf('\n=== STAGE coreg: skipped (DO_COREG_ROIS=0) ===\n');
    else
        if ~haveSessionMean
            sessionMean = resolve_session_mean(funcDir, CFG, false);
        end
        if ~isempty(sessionMean) && (hasRoiL || hasRoiR)
            fprintf('\n=== STAGE coreg: T2w -> %s ===\n', sessionMean);
            roiList = {};
            if hasRoiL, roiList{end+1} = roiL; end %#ok<AGROW>
            if hasRoiR, roiList{end+1} = roiR; end %#ok<AGROW>
            coreg_reslice_rois(sessionMean, t2w, roiList, CFG);
        else
            fprintf('\n=== STAGE coreg: skipped (missing session mean or ROIs) ===\n');
        end
    end
end

% =====================================================================
%  STAGE smooth — Gaussian smoothing of the preprocessed BOLD files
% =====================================================================
if STAGES.smooth
    if CFG.SMOOTH_FWHM > 0
        fprintf('\n=== STAGE smooth: FWHM = %g mm ===\n', CFG.SMOOTH_FWHM);
        smooth_bold_files(funcDir, CFG);
    else
        fprintf('\n=== STAGE smooth: skipped (SMOOTH_FWHM = 0) ===\n');
    end
end

% --------- Save preprocessing provenance log ----------
save_provenance_log(DER, SUB, SES, CFG, cfgFile, taskBoldFiles, boldPicks, STAGES.list);

fprintf('\nDONE: %s / %s  [stages: %s]\n\n', SUB, SES, strjoin(STAGES.list, ','));
end

% ============================== HELPERS ==============================

function out = newest_match(folder, regexPattern)
d = dir(folder);
bestTime = -Inf;
out = '';
for i=1:numel(d)
    if d(i).isdir, continue; end
    if ~isempty(regexp(d(i).name, regexPattern, 'once'))
        if d(i).datenum > bestTime
            bestTime = d(i).datenum;
            out = fullfile(folder, d(i).name);
        end
    end
end
end

% ===================== CHOOSING WHICH RUN TO USE =====================
%
% A session can hold several runs of the same scan. The run numbers are
% counted per modality and carry no cross-modality meaning: if the subject
% climbed out of the scanner between sequences, run-2 of the T2w and run-2 of
% a task are from different visits and must not be paired.
%
% So the BOLD run being preprocessed is the anchor, and every other input is
% matched to it: an explicit choice in run_selection.tsv first, otherwise the
% run acquired closest in time, otherwise a refusal (or a loud guess).

function [picks, anchor] = select_task_bolds(rawFuncDir, SUB, SES, CFG, runSel, required)
% SELECT_TASK_BOLDS  One BOLD run per task label, and the anchor time.
%
% Returns a struct array with one entry per CFG.TASK_LABELS:
%   .label .path .runLabel .acqSec .reason
% and the anchor, which is what the anatomy and fieldmaps are matched against:
%   .secs  every chosen functional run
%   .sec   the earliest of them, and .name the task it belongs to
%
% Tasks are resolved in two passes. A task with only one run needs no decision
% and pins down where the functional block is; those go first, and the tasks
% that do need a decision are then matched against them. Resolving in config
% order instead would leave the first task with nothing to match to — exactly
% the case where a session has two runs of task-run1 and one of task-run2.

nTasks = numel(CFG.TASK_LABELS);
picks = repmat(struct('label', '', 'path', '', 'runLabel', '', ...
                      'acqSec', NaN, 'reason', ''), nTasks, 1);

anchor = struct('sec', NaN, 'secs', NaN, 'name', 'the functional data');

% --------- Gather the candidates for every task ----------
cands = cell(nTasks, 1);
for ii = 1:nTasks
    tLabel = CFG.TASK_LABELS{ii};
    picks(ii).label = tLabel;

    cands{ii} = bids_list_runs(rawFuncDir, struct('sub', SUB, 'ses', SES, ...
                                                  'suffix', 'bold', 'task', tLabel, ...
                                                  'prefix', ''));

    if isempty(cands{ii}) && required
        error('run_spm_preproc:NoBold', [ ...
            'No BOLD file for task-%s in %s\n' ...
            'Expected something like %s_%s_task-%s[_run-N]_bold.nii[.gz]\n' ...
            'Check the TASKS setting in pipeline_config.cfg against the filenames.'], ...
            tLabel, rawFuncDir, SUB, SES, tLabel);
    end
end

% --------- Pass 1: the tasks that decide themselves ----------
% Either only one run exists, or run_selection.tsv already names one.
resolved = false(nTasks, 1);
anchorSecs = [];

for ii = 1:nTasks
    if isempty(cands{ii}), continue; end
    tLabel = picks(ii).label;
    selected = bids_selection_for(runSel, ['task-' tLabel], tLabel);

    if numel(cands{ii}) > 1 && isempty(selected)
        continue;   % needs an anchor; pass 2
    end

    [p, info] = bids_pick_run(cands{ii}, struct( ...
        'label',    sprintf('BOLD task-%s', tLabel), ...
        'selected', selected, ...
        'gapSec',   CFG.RUN_MATCH_GAP_SEC, ...
        'strict',   CFG.STRICT_RUN_MATCHING));

    picks(ii) = store_pick(picks(ii), p, info);
    resolved(ii) = true;
    if ~isnan(info.acqSec), anchorSecs(end+1) = info.acqSec; end %#ok<AGROW>
    fprintf('  BOLD task-%-8s %-10s %s\n', tLabel, info.runLabel, info.reason);
end

% --------- Pass 2: the rest, matched to what pass 1 settled ----------
for ii = 1:nTasks
    if resolved(ii) || isempty(cands{ii}), continue; end
    tLabel = picks(ii).label;

    [p, info] = bids_pick_run(cands{ii}, struct( ...
        'label',      sprintf('BOLD task-%s', tLabel), ...
        'selected',   '', ...
        'anchorSec',  anchor_vector(anchorSecs), ...
        'anchorName', 'the other functional runs', ...
        'gapSec',     CFG.RUN_MATCH_GAP_SEC, ...
        'strict',     CFG.STRICT_RUN_MATCHING));

    picks(ii) = store_pick(picks(ii), p, info);
    if ~isnan(info.acqSec), anchorSecs(end+1) = info.acqSec; end %#ok<AGROW>
    fprintf('  BOLD task-%-8s %-10s %s\n', tLabel, info.runLabel, info.reason);
end

% --------- The anchor everything else is matched against ----------
secs = [picks.acqSec];
secs = secs(~isnan(secs));
if ~isempty(secs)
    anchor.secs = secs;
    anchor.sec  = min(secs);
    idx = 1;
    for ii = 1:nTasks
        if ~isnan(picks(ii).acqSec) && picks(ii).acqSec == anchor.sec
            idx = ii;
            break;
        end
    end
    anchor.name = sprintf('task-%s %s', picks(idx).label, picks(idx).runLabel);

    % Functional runs spread over more than one visit: the anatomy can only
    % match one of them, so say so before anything is coregistered.
    if max(secs) - min(secs) > CFG.RUN_MATCH_GAP_SEC
        warning('run_spm_preproc:TasksSpanBlocks', [ ...
            'The chosen functional runs span %.0f minutes, more than\n' ...
            'RUN_MATCH_GAP_MIN. They were probably not acquired in one visit to the\n' ...
            'scanner, so one session mean covers two head positions and the T2w can\n' ...
            'only match one of them. Check run_selection.tsv for %s / %s.'], ...
            (max(secs) - min(secs))/60, SUB, SES);
    end
end
end

function pick = store_pick(pick, p, info)
pick.path     = p;
pick.runLabel = info.runLabel;
pick.acqSec   = info.acqSec;
pick.reason   = info.reason;
end

function v = anchor_vector(secs)
if isempty(secs)
    v = NaN;
else
    v = secs;
end
end

function t2wFile = select_t2w(rawAnatDir, SUB, SES, CFG, runSel, anchor)
% SELECT_T2W  The structural to coregister to the session mean.

cands = bids_list_runs(rawAnatDir, struct('sub', SUB, 'ses', SES, ...
                                          'suffix', 'T2w', 'prefix', ''));

if isempty(cands)
    error('run_spm_preproc:NoT2w', [ ...
        'No T2w file in %s\n' ...
        'Expected %s_%s[_run-N]_T2w.nii[.gz]'], rawAnatDir, SUB, SES);
end

[t2wFile, info] = bids_pick_run(cands, struct( ...
    'label',      'T2w', ...
    'selected',   bids_selection_for(runSel, 'T2w'), ...
    'anchorSec',  anchor.secs, ...
    'anchorName', anchor.name, ...
    'gapSec',     CFG.RUN_MATCH_GAP_SEC, ...
    'strict',     CFG.STRICT_RUN_MATCHING));

fprintf('  T2w      %-10s %s\n', info.runLabel, info.reason);
end

function [phasediff, magFile] = select_phasediff(rawFmapDir, SUB, SES, CFG, runSel, anchor)
% SELECT_PHASEDIFF  A GRE phasediff plus the magnitude1 acquired with it.
%
% The magnitude has to come from the same acquisition as the phase difference,
% not merely from the same session: it is what the phase is unwrapped against.
% So the phasediff is chosen first and the magnitude is taken from the file
% carrying the same entities; only if there is no such file does it fall back
% to matching by acquisition time.

cands = bids_list_runs(rawFmapDir, struct('sub', SUB, 'ses', SES, ...
                                          'suffix', 'phasediff', 'prefix', ''));
if isempty(cands)
    error('run_spm_preproc:NoPhasediff', [ ...
        'No phasediff in %s\n' ...
        'PREPROC_MODE=realign_unwarp needs %s_%s[_run-N]_phasediff.nii[.gz]\n' ...
        'plus its magnitude1.'], rawFmapDir, SUB, SES);
end

[phasediff, info] = bids_pick_run(cands, struct( ...
    'label',      'phasediff', ...
    'selected',   bids_selection_for(runSel, 'phasediff', 'fieldmap'), ...
    'anchorSec',  anchor.secs, ...
    'anchorName', anchor.name, ...
    'gapSec',     CFG.RUN_MATCH_GAP_SEC, ...
    'strict',     CFG.STRICT_RUN_MATCHING));
fprintf('  phasediff %-9s %s\n', info.runLabel, info.reason);

magFile = companion_file(phasediff, rawFmapDir, 'magnitude1', SUB, SES, CFG, runSel, anchor);
assert(~isempty(magFile), [ ...
    'Found the phasediff but no magnitude1 to unwrap it against.\n' ...
    '  phasediff: %s\n' ...
    '  Looked in: %s\n' ...
    'SPM needs the magnitude acquired with this phase difference. If the\n' ...
    'session has none, build a fieldmap from the T1w instead:\n' ...
    '  bash make_fieldmaps.sh --sub %s --ses %s\n' ...
    'then set PREPROC_MODE=precalc_fieldmap.'], phasediff, rawFmapDir, SUB, SES);
end

function [fieldmapFile, magFile] = select_precalc_fieldmap(rawFmapDir, derivFmapDir, SUB, SES, CFG, runSel, anchor)
% SELECT_PRECALC_FIELDMAP  A ready-made B0 map and the magnitude it belongs to.
%
% make_fieldmaps.sh writes these into the fieldmaps derivative; older datasets
% may still have them next to the raw fmap files, so both are searched.

searchDirs = {derivFmapDir, rawFmapDir};

fmPat  = CFG.FIELDMAP_PATTERN;            % '_fieldmap'
magPat = CFG.FIELDMAP_MAGNITUDE_PATTERN;  % '_magnitude'

cands = bids_list_runs(searchDirs, struct('sub', SUB, 'ses', SES, ...
    'regexp', ['.*' regexptranslate('escape', fmPat) '\.nii(\.gz)?$'], 'prefix', ''));

assert(~isempty(cands), [ ...
    'No fieldmap found for %s / %s.\n' ...
    '  Looked in: %s\n' ...
    '             %s\n' ...
    'Expected a file like %s_%s%s.nii[.gz]\n' ...
    'Create one with:  bash make_fieldmaps.sh --sub %s --ses %s\n' ...
    '(or point FIELDMAP_PATTERN at whatever your fieldmaps are called)'], ...
    SUB, SES, derivFmapDir, rawFmapDir, SUB, SES, fmPat, SUB, SES);

[fieldmapFile, info] = bids_pick_run(cands, struct( ...
    'label',      'fieldmap', ...
    'selected',   bids_selection_for(runSel, 'fieldmap'), ...
    'anchorSec',  anchor.secs, ...
    'anchorName', anchor.name, ...
    'gapSec',     CFG.RUN_MATCH_GAP_SEC, ...
    'strict',     CFG.STRICT_RUN_MATCHING));
fprintf('  fieldmap %-10s %s\n', info.runLabel, info.reason);

% The magnitude must be the one written next to this fieldmap: same entities,
% same directory, same acquisition.
fmDir  = fileparts(fieldmapFile);
fmEnt  = bids_entities(fieldmapFile);
fmStem = regexprep(fmEnt.stem, [regexptranslate('escape', fmPat) '$'], '');

magCands = bids_list_runs(fmDir, struct( ...
    'regexp', ['^' regexptranslate('escape', fmStem) regexptranslate('escape', magPat) '\.nii(\.gz)?$']));

if isempty(magCands)
    % Fall back to any magnitude for this session, matched by time
    magCands = bids_list_runs(searchDirs, struct('sub', SUB, 'ses', SES, ...
        'regexp', ['.*' regexptranslate('escape', magPat) '\.nii(\.gz)?$'], 'prefix', ''));
end

assert(~isempty(magCands), [ ...
    'Found the fieldmap but not its magnitude image.\n' ...
    '  Fieldmap: %s\n' ...
    '  Expected: %s%s.nii[.gz] next to it\n' ...
    'SPM needs a magnitude in the same space as the fieldmap to mask it and to\n' ...
    'match the VDM to the EPI. make_fieldmaps.sh writes one beside every\n' ...
    'fieldmap it creates.'], fieldmapFile, fmStem, magPat);

[magFile, magInfo] = bids_pick_run(magCands, struct( ...
    'label',      'fieldmap magnitude', ...
    'selected',   bids_selection_for(runSel, 'magnitude'), ...
    'anchorSec',  bids_acq_time(fieldmapFile), ...
    'anchorName', 'the fieldmap', ...
    'gapSec',     CFG.RUN_MATCH_GAP_SEC, ...
    'strict',     false));
fprintf('  magnitude %-9s %s\n', magInfo.runLabel, magInfo.reason);
end

function reverseFile = select_reverse_pe(RAW, SUB, SES, CFG, runSel, anchor)
% SELECT_REVERSE_PE  The reverse phase-encode EPI topup estimates the field from.
%
% Which directory to search comes from TOPUP_REVERSE_PE_DIR, and what the file
% is called from TOPUP_REVERSE_PE_PATTERN (a substring, so extra entities are
% tolerated). Getting the run wrong here is the worst case of all: the field
% would be estimated between EPIs from two different head positions.

switch CFG.TOPUP_REVERSE_PE_DIR
    case 'fmap', searchDirs = {RAW.fmap};
    case 'func', searchDirs = {RAW.func};
    case 'auto', searchDirs = {RAW.fmap, RAW.func};
    otherwise
        error('Invalid TOPUP_REVERSE_PE_DIR: %s (must be fmap, func, or auto)', ...
              CFG.TOPUP_REVERSE_PE_DIR);
end

pattern = CFG.TOPUP_REVERSE_PE_PATTERN;
cands = bids_list_runs(searchDirs, struct('sub', SUB, 'ses', SES, ...
                                          'contains', pattern, 'prefix', ''));

if isempty(cands)
    % Same fallbacks as before: the usual names for a reverse-PE EPI
    fallbacks = {'dir-[A-Za-z]+_epi', 'task-reverse.*_bold', 'task-reverse', '_epi'};
    for f = 1:numel(fallbacks)
        cands = bids_list_runs(searchDirs, struct('sub', SUB, 'ses', SES, ...
                                                  'regexp', ['.*' fallbacks{f} '.*\.nii(\.gz)?$'], ...
                                                  'prefix', ''));
        if ~isempty(cands)
            fprintf('  (matched via fallback pattern: %s)\n', fallbacks{f});
            break;
        end
    end
end

assert(~isempty(cands), [ ...
    'Could not find the reverse-PE EPI for %s / %s.\n' ...
    '  Pattern:  %s\n' ...
    '  Searched: %s\n' ...
    '  TOPUP_REVERSE_PE_DIR: %s\n' ...
    'Adjust TOPUP_REVERSE_PE_PATTERN in pipeline_config.cfg.'], ...
    SUB, SES, pattern, strjoin(searchDirs, ', '), CFG.TOPUP_REVERSE_PE_DIR);

[reverseFile, info] = bids_pick_run(cands, struct( ...
    'label',      'reverse-PE EPI', ...
    'selected',   bids_selection_for(runSel, 'reverse', 'task-reverse', 'epi'), ...
    'anchorSec',  anchor.secs, ...
    'anchorName', anchor.name, ...
    'gapSec',     CFG.RUN_MATCH_GAP_SEC, ...
    'strict',     CFG.STRICT_RUN_MATCHING));
fprintf('  reverse  %-10s %s\n', info.runLabel, info.reason);
end

function [roiL, roiR, hasRoiL, hasRoiR] = find_rois(searchDirs, SUB, SES, CFG, runSel, anchor)
% FIND_ROIS  Left and right ROI masks for this session.
%
% searchDirs are tried in order, so a mask imported into the derivatives by
% move_rois.sh wins over one still sitting in the raw anat/ directory.
%
% ROI_PATTERN_LEFT/RIGHT are matched as a suffix, which leaves room for the
% extra entities a multi-run session carries. When several masks match, the
% one whose run entity was selected (or, failing that, whose T2w run was
% chosen) is used — a mask drawn on run-1 of the T2w is meaningless on run-2.

roiL = ''; roiR = ''; hasRoiL = false; hasRoiR = false;

if ~CFG.DO_COREG_ROIS, return; end

[roiL, hasRoiL] = pick_roi(searchDirs, SUB, SES, CFG.ROI_PATTERN_LEFT,  'left',  CFG, runSel, anchor);
[roiR, hasRoiR] = pick_roi(searchDirs, SUB, SES, CFG.ROI_PATTERN_RIGHT, 'right', CFG, runSel, anchor);
end

function [roiFile, found] = pick_roi(searchDirs, SUB, SES, pattern, side, CFG, runSel, anchor)
roiFile = ''; found = false;
if isempty(pattern), return; end

escPat = regexptranslate('escape', pattern);

for d = 1:numel(searchDirs)
    cands = bids_list_runs(searchDirs{d}, struct('sub', SUB, 'ses', SES, ...
        'regexp', ['.*' escPat '(\.gz)?$'], 'prefix', ''));
    if isempty(cands), continue; end

    [roiFile, info] = bids_pick_run(cands, struct( ...
        'label',      sprintf('%s ROI mask', side), ...
        'selected',   bids_selection_for(runSel, ['roi-' side], 'roi', 'T2w'), ...
        'anchorSec',  anchor.secs, ...
        'anchorName', anchor.name, ...
        'gapSec',     CFG.RUN_MATCH_GAP_SEC, ...
        'strict',     false));   % a missing sidecar is normal for a hand-drawn mask
    found = true;
    fprintf('  ROI %-5s %-10s %s\n', side, info.runLabel, info.reason);
    return;
end
end

function out = companion_file(refFile, folder, suffix, SUB, SES, CFG, runSel, anchor)
% COMPANION_FILE  The file acquired with refFile, by entities first, time second.
%
% Used for magnitude/phasediff pairs, where "the same session" is not good
% enough: they have to be the same acquisition.

out = '';
refEnt = bids_entities(refFile);

% Same run entity as the reference — this is the reliable case
sameRun = bids_list_runs(folder, struct('sub', SUB, 'ses', SES, ...
                                        'suffix', suffix, 'prefix', ''));
if isempty(sameRun), return; end

for i = 1:numel(sameRun)
    if strcmp(sameRun(i).ent.run, refEnt.run)
        out = sameRun(i).path;
        return;
    end
end

% No entity match: fall back to whichever was acquired closest to it
[out, info] = bids_pick_run(sameRun, struct( ...
    'label',      suffix, ...
    'selected',   bids_selection_for(runSel, suffix), ...
    'anchorSec',  bids_acq_time(refFile), ...
    'anchorName', refEnt.suffix, ...
    'gapSec',     CFG.RUN_MATCH_GAP_SEC, ...
    'strict',     false, ...
    'required',   false));
if ~isempty(out)
    fprintf('  %-9s %-10s %s\n', suffix, info.runLabel, info.reason);
end
end

function write_run_manifest(DER, SUB, SES, CFG, boldPicks, others, anchor)
% WRITE_RUN_MANIFEST  Record which run of each scan this session used.
%
% Written as <ses>/preproc_run_selection.tsv. Two readers need it:
%   - prepare_gridcat_directory.m, to pair each task with its event table even
%     when the event file's name says nothing about runs;
%   - you, six months later, when a session looks odd and the question is
%     which T2w it was coregistered to.
%
% Rows are merged, not replaced. Running the coreg stage on its own decides
% the T2w and the ROI masks but nothing about the BOLD, and the record of
% which BOLD run the realignment used has to survive that.

outFile = fullfile(DER.sesDir, 'preproc_run_selection.tsv');

rows = read_existing_manifest(outFile);

for i = 1:numel(boldPicks)
    if isempty(boldPicks(i).path), continue; end
    key = ['task-' boldPicks(i).label];
    rows.(matlab.lang.makeValidName(key)) = sprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s', ...
        SUB, SES, key, boldPicks(i).runLabel, clock_of(boldPicks(i).acqSec), ...
        boldPicks(i).path, boldPicks(i).reason);
end

types = fieldnames(others);
for i = 1:numel(types)
    f = others.(types{i});
    if isempty(f), continue; end
    ent = bids_entities(f);
    runLabel = ent.run;
    if isempty(runLabel), runLabel = 'no-run'; end
    key = strrep(types{i}, '_', '-');
    rows.(matlab.lang.makeValidName(key)) = sprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s', ...
        SUB, SES, key, runLabel, clock_of(bids_acq_time(f)), f, ...
        'matched to the functional runs');
end

fid = fopen(outFile, 'w');
if fid < 0
    warning('run_spm_preproc:ManifestFailed', 'Could not write %s', outFile);
    return;
end
cleanupFid = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'subject\tsession\ttype\tselected_run\tacq_time\tsource_file\treason\n');
keys = sort(fieldnames(rows));
for i = 1:numel(keys)
    fprintf(fid, '%s\n', rows.(keys{i}));
end
fprintf(fid, '#\tanchor\t%s\t%s\tgap=%gmin\tstrict=%d\t\n', ...
    anchor.name, clock_of(anchor.sec), CFG.RUN_MATCH_GAP_SEC/60, CFG.STRICT_RUN_MATCHING);

fprintf('Run selection written to: %s\n', outFile);
end

function rows = read_existing_manifest(outFile)
% READ_EXISTING_MANIFEST  Previous rows, keyed by type, so they can be kept.
rows = struct();
if ~isfile(outFile), return; end

fid = fopen(outFile, 'r');
if fid == -1, return; end
cleanupFid = onCleanup(@() fclose(fid)); %#ok<NASGU>

lineNo = 0;
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    lineNo = lineNo + 1;
    if lineNo == 1, continue; end                    % header
    if isempty(strtrim(line)), continue; end
    if line(1) == '#', continue; end                 % the anchor footer

    parts = strsplit(line, sprintf('\t'));
    if numel(parts) < 3, continue; end
    key = matlab.lang.makeValidName(strtrim(parts{3}));
    rows.(key) = line;
end
end

function s = clock_of(sec)
if isnan(sec)
    s = 'unknown';
else
    s = sprintf('%02d:%02d:%02d', floor(sec/3600), floor(mod(sec,3600)/60), floor(mod(sec,60)));
end
end

function doReset = deriv_reset_wanted(CFG, STAGES)
% DERIV_RESET_WANTED  Should this run start from an empty derivatives folder?
%
% Only a run that redoes the preprocessing from its beginning may clear the
% previous output. A resumed run (realign,coreg,smooth and friends) reads that
% output back, so wiping it would delete its own input.

switch CFG.DERIV_RESET
    case {'never', 'false', 'no', '0'}
        doReset = false;
        return;
    case {'always', 'true', 'yes', '1'}
        doReset = true;
        return;
end

% 'auto': the first stage this mode actually performs is the from-scratch mark
switch CFG.PREPROC_MODE
    case 'topup',                          entry = 'topup';
    case {'realign_unwarp','precalc_fieldmap'}, entry = 'vdm';
    otherwise,                             entry = 'realign';
end
doReset = STAGES.(entry);
end

function v = cfg_flag(cfgRaw, key, default)
% CFG_FLAG  Read a true/false setting out of the config.
v = default;
if ~isfield(cfgRaw, key), return; end
raw = cfgRaw.(key);
if isempty(raw), return; end
if isnumeric(raw)
    v = raw ~= 0;
    return;
end
v = ismember(lower(strtrim(char(raw))), {'true', 'yes', '1', 'on'});
end

function out = tern(cond, a, b)
if cond, out = a; else, out = b; end
end



function vols = expand_4d(niiFile)
V = spm_vol(niiFile);
vols = cell(numel(V),1);
for k=1:numel(V)
    vols{k} = sprintf('%s,%d', niiFile, k);
end
end

%

function vdm_out = calc_vdm_for_run_to_func(phasemap, magnitude1, epiRef, funcDir, fmapDir, CFG, taskLabel)
% Calculates a VDM matched to epiRef and ensures output ends up in funcDir
% as a deterministic name: vdm_task-<label>.nii
%
% Robust strategy:
% - Don't trust vdmflags.prefix (SPM/FieldMap sometimes ignores or overrides it)
% - Instead: snapshot existing vdm* files (func+fmap), run job, then diff to find newly created VDM.

% ---- Snapshot existing VDM candidates (before) ----
vdmSearchDirs = {funcDir, fmapDir, pwd};
before = list_vdm_candidates(vdmSearchDirs);

% ---- Build FieldMap batch ----
matlabbatch = {};

matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase     = {phasemap};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = {magnitude1};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session.epi = {epiRef};

matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et        = [CFG.TE_SHORT CFG.TE_LONG];
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = CFG.MASKBRAIN;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir   = CFG.BLIPDIR;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert      = CFG.TOTAL_READOUT;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm     = CFG.EPI_BASED_FIELDMAP;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm       = 0;

matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = CFG.VDM_UFLAGS_METHOD;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm   = CFG.VDM_UFLAGS_FWHM;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad    = CFG.VDM_UFLAGS_PAD;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws     = CFG.VDM_UFLAGS_WS;

matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {fullfile(spm('Dir'),'toolbox','FieldMap','T1.nii')};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;

matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm      = 1;  % write VDM matched to epiRef
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname      = sprintf('task-%s', taskLabel);
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;

% These avoid "incomplete module inputs" errors on some installs
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat      = {''};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;

% vdmflags: keep them, but DON'T depend on prefix to identify output
% (prefix may get ignored/overridden, but flags still help define output behavior)
vdmflags = CFG.VDMFLAGS_BASE;
% Leave prefix at whatever you set in CFG.VDMFLAGS_BASE (e.g., 'vdm_') or default-like 'vdm5_'
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.vdmflags = vdmflags;

% ---- Run job (run from funcDir to make outputs more likely to land there) ----
oldpwd = pwd;
cleanup = onCleanup(@() cd(oldpwd));
cd(funcDir);

spm_jobman('run', matlabbatch);

% ---- Snapshot after, diff to find NEW VDM(s) ----
after = list_vdm_candidates([vdmSearchDirs, {pwd}]);
newFiles = diff_vdm_candidates(before, after);

assert(~isempty(newFiles), 'VDM not found after FieldMap run for task-%s. Check where FieldMap writes outputs.', taskLabel);

% If multiple candidates were created, take the newest one
newest = pick_newest_file(newFiles);

% ---- Move/rename deterministically into funcDir ----
vdm_out = fullfile(funcDir, sprintf('vdm_task-%s.nii', taskLabel));

% If output is .img/.hdr pair, convert handling:
[~,~,ext] = fileparts(newest);
if strcmpi(ext,'.img')
    % Move both .img and .hdr, but keep .img name; SPM can read Analyze too.
    target_img = fullfile(funcDir, sprintf('vdm_task-%s.img', taskLabel));
    target_hdr = fullfile(funcDir, sprintf('vdm_task-%s.hdr', taskLabel));
    if isfile(target_img), delete(target_img); end
    if isfile(target_hdr), delete(target_hdr); end
    movefile(newest, target_img);
    hdr = strrep(newest,'.img','.hdr');
    if isfile(hdr), movefile(hdr, target_hdr); end
    % Return .img path in this case
    vdm_out = target_img;
else
    % Assume NIfTI .nii
    if isfile(vdm_out), delete(vdm_out); end
    movefile(newest, vdm_out);
end

fprintf('Renamed VDM -> %s\n', vdm_out);

end

% ---------------- helper: list VDM candidates ----------------
function files = list_vdm_candidates(folders)
% LIST_VDM_CANDIDATES  Snapshot the vdm* files across one or more folders.
%
% Returns a struct array with .path, .datenum and .bytes. Recording the
% timestamp and size (not just the path) means a VDM that SPM overwrote in
% place still registers as new output — which happens when a stale vdm5_* from
% an earlier manual run is already sitting in fmap/.

if ischar(folders), folders = {folders}; end
files = struct('path', {}, 'datenum', {}, 'bytes', {});
seen = {};

for f = 1:numel(folders)
    folder = folders{f};
    if isempty(folder) || ~isfolder(folder), continue; end
    d = dir(folder);
    for i = 1:numel(d)
        if d(i).isdir, continue; end
        n = d(i).name;

        % FieldMap typically uses vdm*.nii or vdm*.img/.hdr
        if isempty(regexp(n,'^vdm.*\.(nii|img)$','once')), continue; end

        p = fullfile(folder, n);
        if ismember(p, seen), continue; end
        seen{end+1} = p; %#ok<AGROW>
        files(end+1) = struct('path', p, 'datenum', d(i).datenum, 'bytes', d(i).bytes); %#ok<AGROW>
    end
end
end

% ---------------- helper: VDM files that appeared or changed ----------------
function newFiles = diff_vdm_candidates(before, after)
% DIFF_VDM_CANDIDATES  Paths in `after` that are new, or were rewritten.

newFiles = {};
beforePaths = {before.path};

for i = 1:numel(after)
    idx = find(strcmp(after(i).path, beforePaths), 1);
    if isempty(idx)
        newFiles{end+1,1} = after(i).path; %#ok<AGROW>
    elseif after(i).datenum > before(idx).datenum || after(i).bytes ~= before(idx).bytes
        newFiles{end+1,1} = after(i).path; %#ok<AGROW>
    end
end
end

% ---------------- helper: pick newest file ----------------
function f = pick_newest_file(fileList)
assert(~isempty(fileList));
bestT = -Inf;
f = fileList{1};

for i = 1:numel(fileList)
    info = dir(fileList{i});
    if ~isempty(info) && info.datenum > bestT
        bestT = info.datenum;
        f = fileList{i};
    end
end
end

function run_realign_unwarp_multi(taskVols, taskVDM, CFG)
% ONE Realign&Unwarp batch, multiple data blocks.
% taskVols: cell array, each element is cellstr of scans for that run
% taskVDM : cell array, each element is full path to vdm for that run
% CFG     : configuration struct with realign/unwarp parameters

n = numel(taskVols);
assert(numel(taskVDM) == n, 'taskVols and taskVDM must have same length.');

matlabbatch = {};
for i = 1:n
    matlabbatch{1}.spm.spatial.realignunwarp.data(i).scans  = taskVols{i};
    matlabbatch{1}.spm.spatial.realignunwarp.data(i).pmscan = {taskVDM{i}};
end

% Estimation options (from config)
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = CFG.REALIGN_QUALITY;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep     = CFG.REALIGN_SEP;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm    = CFG.REALIGN_FWHM;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm     = CFG.REALIGN_RTM;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = CFG.REALIGN_EINTERP;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap   = CFG.REALIGN_EWRAP;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight  = '';

% Unwarp options (from config)
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = CFG.UNWARP_BASFCN;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = CFG.UNWARP_REGORDER;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = CFG.UNWARP_LAMBDA;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = CFG.UNWARP_JM;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = CFG.UNWARP_FOT;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = CFG.UNWARP_SOT;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = CFG.UNWARP_UWFWHM;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = CFG.UNWARP_REM;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = CFG.UNWARP_NOI;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = CFG.UNWARP_EXPROUND;

% Reslice options (from config)
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = CFG.WHICHUNWARP;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = CFG.RESLICE_INTERP;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = CFG.RESLICE_WRAP;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = CFG.RESLICE_MASK;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = CFG.RESLICE_PREFIX;

spm_jobman('run', matlabbatch);
end

function run_realign_only_multi(taskVols, CFG)
% Realign only (no unwarp/VDM).
% taskVols: cell array, each element is cellstr of scans for that run
% CFG     : configuration struct with realign parameters

n = numel(taskVols);

matlabbatch = {};

% Data specification
for i = 1:n
    matlabbatch{1}.spm.spatial.realign.estwrite.data{i} = taskVols{i};
end

% Estimation options (from config)
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = CFG.REALIGN_QUALITY;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = CFG.REALIGN_SEP;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = CFG.REALIGN_FWHM;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = CFG.REALIGN_RTM;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = CFG.REALIGN_EINTERP;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = CFG.REALIGN_EWRAP;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';

% Reslice options
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];  % Write resliced + mean
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = CFG.RESLICE_INTERP;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = CFG.RESLICE_WRAP;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = CFG.RESLICE_MASK;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = CFG.RESLICE_PREFIX;

spm_jobman('run', matlabbatch);
end

function make_session_mean(meanu_imgs, outFile)
% Average the run-level mean images only (fast)
matlabbatch = {};
matlabbatch{1}.spm.util.imcalc.input = meanu_imgs;
matlabbatch{1}.spm.util.imcalc.output = spm_file(outFile,'filename'); % just filename
matlabbatch{1}.spm.util.imcalc.outdir = {spm_file(outFile,'path')};

expr_terms = arrayfun(@(k)sprintf('i%d',k), 1:numel(meanu_imgs), 'UniformOutput', false);
matlabbatch{1}.spm.util.imcalc.expression = sprintf('(%s)/%d', strjoin(expr_terms,'+'), numel(meanu_imgs));

matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 16;

spm_jobman('run', matlabbatch);
end

function coreg_reslice_rois(refEPI, t2w, roiFiles, CFG)
% Coregister: source = T2w, ref = mean EPI
% Reslice ROIs using the estimated transform (NN interp for masks).
% CFG is required for coregistration parameters
matlabbatch = {};

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {refEPI};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {t2w};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = roiFiles(:);
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = CFG.COREG_COST_FUN;
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = CFG.COREG_SEP;
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = ...
    [0.0200 0.0200 0.0200 0.0010 0.0010 0.0010 0.0100 0.0100 0.0100 0.0010 0.0010 0.0010];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = CFG.COREG_FWHM;

% Reslicing options: NN for masks
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = CFG.COREG_ROI_INTERP; % 0 = nearest neighbor for masks
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

spm_jobman('run', matlabbatch);
end


function smooth_bold_files(funcDir, CFG)
% SMOOTH_BOLD_FILES  Apply Gaussian smoothing to preprocessed BOLD files.
%
% Finds all resliced BOLD files (u*_bold*.nii) in funcDir, excluding mean
% images, and smooths them using SPM's smooth module.
% Output files get SMOOTH_PREFIX prepended (e.g. u -> su).

prefix = CFG.RESLICE_PREFIX;  % typically 'u'

% Find preprocessed BOLD files but NOT mean images
d = dir(funcDir);
boldFiles = {};
for i = 1:numel(d)
    if d(i).isdir, continue; end
    % Match files starting with the reslice prefix, containing _bold, ending .nii
    % but exclude mean images (meanu*, mean*)
    if ~isempty(regexp(d(i).name, ['^' prefix '.*_bold.*\.nii$'], 'once')) && ...
       isempty(regexp(d(i).name, '^mean', 'once'))
        boldFiles{end+1} = fullfile(funcDir, d(i).name); %#ok<AGROW>
    end
end

if isempty(boldFiles)
    warning('smooth_bold_files:NoFiles', 'No preprocessed BOLD files found to smooth in %s', funcDir);
    return;
end

fprintf('  Found %d BOLD files to smooth.\n', numel(boldFiles));

% Ensure isotropic kernel
fwhm = CFG.SMOOTH_FWHM;
if isscalar(fwhm)
    fwhm = [fwhm fwhm fwhm];
end

% Build and run SPM smooth batch
matlabbatch = {};
matlabbatch{1}.spm.spatial.smooth.data   = boldFiles(:);
matlabbatch{1}.spm.spatial.smooth.fwhm   = fwhm;
matlabbatch{1}.spm.spatial.smooth.dtype  = 0;   % same datatype as input
matlabbatch{1}.spm.spatial.smooth.im     = 0;   % no implicit masking
matlabbatch{1}.spm.spatial.smooth.prefix = CFG.SMOOTH_PREFIX;

spm_jobman('run', matlabbatch);
fprintf('  Smoothing complete. Output prefix: %s\n', CFG.SMOOTH_PREFIX);
end

% ================= PRECALCULATED FIELDMAP HELPERS =================

function hzFile = fieldmap_to_hz(fieldmapFile, funcDir, CFG)
% FIELDMAP_TO_HZ  Write a copy of the fieldmap in Hz for SPM's FieldMap toolbox.
%
% fsl_prepare_fieldmap writes rad/s; SPM's "Precalculated FieldMap" input is
% documented as Hz, so rad/s is divided by 2*pi. A fieldmap already in Hz is
% copied through unchanged, so the rest of the code has one path to follow.
%
% Doing the unit conversion here — and letting SPM do the Hz -> voxel-shift
% conversion with its own internal scaling (from TOTAL_READOUT_MS and
% BLIP_DIRECTION) — keeps this code out of the business of guessing what
% scaling SPM's VDM files use.

V = spm_vol(fieldmapFile);
if numel(V) > 1, V = V(1); end
dat = spm_read_vols(V);

switch CFG.FIELDMAP_UNITS
    case {'rad/s','rads','rad_per_s'}
        dat = dat / (2*pi);
        fprintf('  Converted fieldmap rad/s -> Hz (divided by 2*pi).\n');
    case 'hz'
        fprintf('  Fieldmap is already in Hz; no conversion.\n');
    otherwise
        error('Unknown FIELDMAP_UNITS: %s (must be rad/s or Hz)', CFG.FIELDMAP_UNITS);
end

fprintf('  Fieldmap range: %.1f to %.1f Hz\n', min(dat(:)), max(dat(:)));

hzFile = fullfile(funcDir, 'fieldmap_hz.nii');
Vout = V;
Vout.fname = hzFile;
Vout.dt = [spm_type('float32') 0];
Vout.pinfo = [1;0;0];
Vout.descrip = 'B0 fieldmap in Hz (for SPM FieldMap)';
spm_write_vol(Vout, dat);
end

function vdm_out = calc_vdm_from_fieldmap(fieldmapHz, magFile, epiRef, funcDir, fmapDir, CFG, taskLabel)
% CALC_VDM_FROM_FIELDMAP  Build a VDM from a precalculated fieldmap (in Hz).
%
% Uses SPM's FieldMap toolbox "Precalculated FieldMap" input, so SPM performs
% the fieldmap -> voxel displacement conversion itself from TOTAL_READOUT_MS
% (tert) and BLIP_DIRECTION, and matches the VDM to this run's EPI geometry.
%
% Same output contract as calc_vdm_for_run_to_func: the VDM ends up in funcDir
% as vdm_task-<label>.nii, found by diffing the vdm* files before and after.

% ---- Snapshot existing VDM candidates (before) ----
vdmSearchDirs = {funcDir, fmapDir, pwd};
before = list_vdm_candidates(vdmSearchDirs);

% ---- Build FieldMap batch ----
matlabbatch = {};

matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.precalcfieldmap.precalcfieldmap = {fieldmapHz};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.precalcfieldmap.magfieldmap     = {magFile};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session.epi = {epiRef};

% Echo times are unused for a precalculated fieldmap (no phase to unwrap), but
% the batch still requires the defaults structure to be complete.
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et        = [CFG.TE_SHORT CFG.TE_LONG];
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = CFG.MASKBRAIN;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir   = CFG.BLIPDIR;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert      = CFG.TOTAL_READOUT;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm     = CFG.EPI_BASED_FIELDMAP;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm       = 0;

matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = CFG.VDM_UFLAGS_METHOD;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm   = CFG.VDM_UFLAGS_FWHM;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad    = CFG.VDM_UFLAGS_PAD;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws     = CFG.VDM_UFLAGS_WS;

matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {fullfile(spm('Dir'),'toolbox','FieldMap','T1.nii')};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;

matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm      = 1;  % write VDM matched to epiRef
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname      = sprintf('task-%s', taskLabel);
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat          = {''};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat     = 0;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.vdmflags      = CFG.VDMFLAGS_BASE;

% ---- Run job (from funcDir, so outputs are more likely to land there) ----
oldpwd = pwd;
cleanup = onCleanup(@() cd(oldpwd)); %#ok<NASGU>
cd(funcDir);

spm_jobman('run', matlabbatch);

% ---- Snapshot after, diff to find NEW (or rewritten) VDM(s) ----
after = list_vdm_candidates([vdmSearchDirs, {pwd}]);
newFiles = diff_vdm_candidates(before, after);

assert(~isempty(newFiles), [ ...
    'VDM not found after FieldMap run for task-%s.\n' ...
    'Check where SPM FieldMap wrote its output (looked in %s, %s).'], ...
    taskLabel, funcDir, fmapDir);

newest = pick_newest_file(newFiles);

% ---- Move/rename deterministically into funcDir ----
vdm_out = fullfile(funcDir, sprintf('vdm_task-%s.nii', taskLabel));

[~,~,ext] = fileparts(newest);
if strcmpi(ext,'.img')
    target_img = fullfile(funcDir, sprintf('vdm_task-%s.img', taskLabel));
    target_hdr = fullfile(funcDir, sprintf('vdm_task-%s.hdr', taskLabel));
    if isfile(target_img), delete(target_img); end
    if isfile(target_hdr), delete(target_hdr); end
    movefile(newest, target_img);
    hdr = strrep(newest,'.img','.hdr');
    if isfile(hdr), movefile(hdr, target_hdr); end
    vdm_out = target_img;
else
    if isfile(vdm_out), delete(vdm_out); end
    movefile(newest, vdm_out);
end

fprintf('Renamed VDM -> %s\n', vdm_out);
end

% ===================== FSL TOPUP HELPER FUNCTIONS =====================

function run_fsl_topup(boldFile, reversePE_epi, topupPrefix, CFG)
% RUN_FSL_TOPUP  Estimate distortion field using FSL topup.
%
% Creates an acqparams.txt file, merges the first BOLD volume with the
% reverse-PE EPI, and runs FSL topup to estimate the fieldmap.
%
% Inputs:
%   boldFile      - Path to the first task's 4D BOLD (uses volume 1)
%   reversePE_epi - Path to the reverse phase-encode EPI in fmap/
%   topupPrefix   - Output prefix for topup results (e.g., funcDir/topup_results)
%   CFG           - Config struct with TOPUP_* fields

funcDir = fileparts(topupPrefix);

% --- Step 1: Extract first volume from BOLD as the "forward" EPI ---
forwardEPI = fullfile(funcDir, 'topup_forward_b0.nii.gz');
cmd = sprintf('fslroi %s %s 0 1', boldFile, forwardEPI);
fprintf('  Extracting first BOLD volume: %s\n', cmd);
run_shell(cmd);

% Handle reverse-PE EPI: might be 4D, take first volume
reverseInfo = nifti_info_via_fsl(reversePE_epi);
if reverseInfo.nvols > 1
    reverseB0 = fullfile(funcDir, 'topup_reverse_b0.nii.gz');
    cmd = sprintf('fslroi %s %s 0 1', reversePE_epi, reverseB0);
    fprintf('  Extracting first reverse-PE volume: %s\n', cmd);
    run_shell(cmd);
else
    reverseB0 = reversePE_epi;
end

% --- Step 2: Merge forward + reverse into a single 4D file ---
mergedFile = fullfile(funcDir, 'topup_merged_b0.nii.gz');
cmd = sprintf('fslmerge -t %s %s %s', mergedFile, forwardEPI, reverseB0);
fprintf('  Merging forward + reverse EPIs: %s\n', cmd);
run_shell(cmd);

% --- Step 3: Create acqparams.txt ---
acqparamsFile = fullfile(funcDir, 'topup_acqparams.txt');
write_acqparams(acqparamsFile, CFG);
fprintf('  Acqparams file: %s\n', acqparamsFile);

% --- Step 4: Run FSL topup ---
topupCmd = sprintf('topup --imain=%s --datain=%s --out=%s --fout=%s_field --iout=%s_corrected', ...
    mergedFile, acqparamsFile, topupPrefix, topupPrefix, topupPrefix);

if ~isempty(CFG.TOPUP_CONFIG)
    topupCmd = sprintf('%s --config=%s', topupCmd, CFG.TOPUP_CONFIG);
end

% Add verbose flag for logging
topupCmd = sprintf('%s --verbose', topupCmd);

fprintf('  Running FSL topup:\n    %s\n', topupCmd);
run_shell(topupCmd);
fprintf('  Topup estimation complete.\n');

% Clean up temporary files
if isfile(forwardEPI), delete(forwardEPI); end
if ~strcmp(reverseB0, reversePE_epi) && isfile(reverseB0), delete(reverseB0); end
if isfile(mergedFile), delete(mergedFile); end

end

function outputFile = run_fsl_applytopup(boldFile, topupPrefix, acqparamsFile, outputFile, imainIndex, CFG)
% RUN_FSL_APPLYTOPUP  Apply topup distortion correction to a 4D BOLD file.
%
% Inputs:
%   boldFile      - Input 4D BOLD NIfTI
%   topupPrefix   - Prefix from run_fsl_topup (contains fieldcoefs, movpar)
%   acqparamsFile - acqparams.txt topup was run with (see resolve_acqparams)
%   outputFile    - Desired output path (must end in .nii)
%   imainIndex    - Index into acqparams.txt for this image (1 = forward PE)
%   CFG           - Config struct with TOPUP_INTERP
%
% Returns:
%   outputFile  - Actual path to the corrected file (.nii, decompressed)
%
% Notes:
%   - FSL --out expects a basename (no .nii extension). We strip it.
%   - FSL always writes .nii.gz output. We decompress to .nii for SPM.
%   - We use --method=jac (Jacobian modulation) for single-image correction.
%     The alternative --method=lsr (least-squares) requires BOTH forward and
%     reverse PE images passed together, which doesn't apply here.

assert(isfile(acqparamsFile), 'acqparams file not found: %s', acqparamsFile);

% FSL --out expects a basename without extension
outBase = regexprep(outputFile, '\.nii(\.gz)?$', '');

cmd = sprintf('applytopup --imain=%s --topup=%s --datain=%s --inindex=%d --method=jac --interp=%s --out=%s', ...
    boldFile, topupPrefix, acqparamsFile, imainIndex, CFG.TOPUP_INTERP, outBase);

fprintf('  Applying topup to %s:\n    %s\n', boldFile, cmd);
run_shell(cmd);

% FSL writes <basename>.nii.gz — decompress to .nii for SPM
outGz = [outBase '.nii.gz'];
outputFile = [outBase '.nii'];

if isfile(outGz)
    fprintf('  Decompressing applytopup output for SPM...\n');
    gunzip(outGz);
    delete(outGz);
end

assert(isfile(outputFile), 'applytopup output not found: %s (checked %s too)', outputFile, outGz);
fprintf('  applytopup complete -> %s\n', outputFile);
end

function vdm_out = convert_topup_to_vdm(topupPrefix, epiRef, funcDir, CFG, taskLabel)
% CONVERT_TOPUP_TO_VDM  Convert FSL topup field to SPM-compatible VDM.
%
% The topup field (in Hz) is converted to a voxel displacement map (in mm)
% using the formula: VDM = fieldmap_Hz * totalReadoutTime * voxelSize_PE
%
% The VDM is then written as a NIfTI in funcDir for use with SPM Realign&Unwarp.

fieldFile = sprintf('%s_field.nii.gz', topupPrefix);
if ~isfile(fieldFile)
    % Topup output made outside this pipeline may not be compressed
    fieldFile = sprintf('%s_field.nii', topupPrefix);
end
assert(isfile(fieldFile), 'Topup field not found: %s_field.nii[.gz]', topupPrefix);

% Convert to .nii for SPM (FSL defaults to .nii.gz)
fieldNii = fullfile(funcDir, sprintf('topup_field_task-%s.nii', taskLabel));
fieldNiiGz = [fieldNii '.gz'];
cmd = sprintf('fslmaths %s %s -odt float', fieldFile, fieldNiiGz);
run_shell(cmd);
% Decompress .nii.gz -> .nii
if isfile(fieldNiiGz)
    gunzip(fieldNiiGz);
    delete(fieldNiiGz);
end
assert(isfile(fieldNii), 'fslmaths output not found: %s', fieldNii);

% Read the field with SPM
Vfield = spm_vol(fieldNii);
fieldData = spm_read_vols(Vfield);

% Read EPI reference to get geometry (first volume only)
% epiRef may contain ',1' frame specifier; if not, spm_vol on 4D returns
% an array of structs — take only the first element.
epiPath = regexp(epiRef, '^[^,]+', 'match', 'once');
Vepi = spm_vol(epiPath);
if numel(Vepi) > 1
    Vepi = Vepi(1);
end

% Determine PE direction index and voxel size along that axis
peDir = CFG.TOPUP_PE_DIR_BOLD;
switch peDir
    case {'x', 'x-'}
        peDim = 1;
    case {'y', 'y-'}
        peDim = 2;
    case {'z', 'z-'}
        peDim = 3;
    otherwise
        error('Unknown PE direction: %s', peDir);
end
voxSize_PE = abs(Vepi.mat(peDim, peDim));

% Convert Hz fieldmap to VDM (voxel shifts in the PE direction)
% VDM = field_Hz * totalReadoutTime_s (gives shift in voxels),
% then multiply by voxel size to get mm displacement
readout_sec = CFG.TOPUP_READOUT_SEC;
vdmData = fieldData * readout_sec * voxSize_PE;

% If PE direction is negative, flip the sign
if endsWith(peDir, '-')
    vdmData = -vdmData;
end

% Write VDM using field header (same geometry)
vdm_out = fullfile(funcDir, sprintf('vdm_task-%s.nii', taskLabel));
Vout = Vfield;
Vout.fname = vdm_out;
Vout.dt = [spm_type('float32') 0];
Vout.descrip = 'VDM from FSL topup';
spm_write_vol(Vout, vdmData);

% Clean up intermediate field file
if isfile(fieldNii), delete(fieldNii); end

end

function vec = pe_dir_to_vector(peDir)
% PE_DIR_TO_VECTOR  Convert FSL PE direction string to acqparams vector.
%   'x'  -> [1 0 0],  'x-' -> [-1 0 0]
%   'y'  -> [0 1 0],  'y-' -> [0 -1 0]
%   'z'  -> [0 0 1],  'z-' -> [0 0 -1]

switch peDir
    case 'x',  vec = [1  0  0];
    case 'x-', vec = [-1 0  0];
    case 'y',  vec = [0  1  0];
    case 'y-', vec = [0 -1  0];
    case 'z',  vec = [0  0  1];
    case 'z-', vec = [0  0 -1];
    otherwise, error('Invalid PE direction: %s (use x, x-, y, y-, z, or z-)', peDir);
end
end

function info = nifti_info_via_fsl(niiFile)
% NIFTI_INFO_VIA_FSL  Get basic NIfTI info using fslinfo.
%   Returns struct with field 'nvols' (number of volumes).

[status, result] = system(sprintf('fslnvols %s', niiFile));
assert(status == 0, 'fslnvols failed for %s: %s', niiFile, result);
info.nvols = str2double(strtrim(result));
end

function run_shell(cmd)
% RUN_SHELL  Execute a shell command and assert success.
[status, result] = system(cmd);
if status ~= 0
    error('Shell command failed (exit %d):\n  Command: %s\n  Output: %s', status, cmd, result);
end
end

% ===================== BIDS JSON VALIDATION =====================

function validate_bids_params(boldFiles, reverseFile, phasediffFile, fieldmapFile, CFG)
% VALIDATE_BIDS_PARAMS  Compare the config against the sidecars of the files
% this run actually chose.
%
% The files are the staged copies in the derivatives, and their JSON sidecars
% came along with them. Checking those rather than "whichever sidecar turns up
% first in the folder" matters as soon as a session has several runs: the run
% being preprocessed is the one whose acquisition parameters have to match the
% config, and a mismatch used to hide behind another run's sidecar.
%
% Nothing here stops the pipeline. The config values are always the ones used,
% because they are sometimes deliberate overrides.

fprintf('\n--- Validating config against the BIDS JSON sidecars ---\n');
nWarnings = 0;

% --- The BOLD run being preprocessed ---
boldJson = '';
for i = 1:numel(boldFiles)
    if isempty(boldFiles{i}), continue; end
    cand = regexprep(boldFiles{i}, '\.nii(\.gz)?$', '.json');
    if isfile(cand)
        boldJson = cand;
        break;
    end
end

bidsPE = '';
bidsPE_hasSign = false;

if ~isempty(boldJson)
    bj = read_json_file(boldJson);
    fprintf('  BOLD JSON: %s\n', boldJson);

    % -- TR --
    if isfield(bj, 'RepetitionTime')
        bidsVal = bj.RepetitionTime;
        if isfield(CFG, 'TR') && abs(CFG.TR - bidsVal) > 0.001
            fprintf('  WARNING: TR mismatch — config: %.4f s, BIDS JSON: %.4f s\n', CFG.TR, bidsVal);
            nWarnings = nWarnings + 1;
        else
            fprintf('  TR: %.4f s (matches)\n', bidsVal);
        end
    end

    % -- Phase encoding direction --
    % BIDS has two PE fields:
    %   PhaseEncodingDirection = full direction with sign (i, i-, j, j-, k, k-)
    %   PhaseEncodingAxis      = axis only, no sign (i, j, k)
    % We use Direction if available, otherwise fall back to Axis.
    if isfield(bj, 'PhaseEncodingDirection') && ~isempty(bj.PhaseEncodingDirection)
        bidsPE = bj.PhaseEncodingDirection;
        bidsPE_hasSign = true;
    elseif isfield(bj, 'PhaseEncodingAxis') && ~isempty(bj.PhaseEncodingAxis)
        bidsPE = bj.PhaseEncodingAxis;
        bidsPE_hasSign = false;
    end
    if ~isempty(bidsPE)
        fslPE = bids_pe_to_fsl(bidsPE);
        bidsAxis = regexprep(fslPE, '-$', '');

        if strcmp(CFG.PREPROC_MODE, 'topup') && ~isempty(fslPE)
            configAxis = regexprep(CFG.TOPUP_PE_DIR_BOLD, '-$', '');  % 'y-' -> 'y'
            if bidsPE_hasSign
                if ~strcmp(fslPE, CFG.TOPUP_PE_DIR_BOLD)
                    fprintf('  WARNING: PE direction mismatch — config TOPUP_PE_DIR_BOLD: %s, BIDS JSON: %s (=%s in FSL)\n', ...
                        CFG.TOPUP_PE_DIR_BOLD, bidsPE, fslPE);
                    nWarnings = nWarnings + 1;
                else
                    fprintf('  PE direction (BOLD): %s = %s in FSL (matches)\n', bidsPE, fslPE);
                end
            else
                if ~strcmp(bidsAxis, configAxis)
                    fprintf('  WARNING: PE axis mismatch — config TOPUP_PE_DIR_BOLD axis: %s, BIDS JSON: %s (=%s in FSL)\n', ...
                        configAxis, bidsPE, bidsAxis);
                    nWarnings = nWarnings + 1;
                else
                    fprintf('  PE axis (BOLD): %s (matches). Note: BIDS JSON has PhaseEncodingAxis only (no sign) — cannot verify +/- direction.\n', bidsAxis);
                end
            end
        end

        if strcmp(CFG.PREPROC_MODE, 'realign_unwarp') && ~isempty(fslPE)
            if bidsPE_hasSign
                % Check blip direction consistency: j -> +1, j- -> -1
                expectedBlip = 1;
                if endsWith(bidsPE, '-'), expectedBlip = -1; end
                if CFG.BLIPDIR ~= expectedBlip
                    fprintf('  WARNING: BLIP_DIRECTION mismatch — config: %d, BIDS JSON PE=%s suggests: %d\n', ...
                        CFG.BLIPDIR, bidsPE, expectedBlip);
                    nWarnings = nWarnings + 1;
                else
                    fprintf('  BLIP_DIRECTION: %d (consistent with BIDS PE=%s)\n', CFG.BLIPDIR, bidsPE);
                end
            else
                fprintf('  BIDS JSON has PhaseEncodingAxis=%s only (no sign) — cannot verify BLIP_DIRECTION.\n', bidsPE);
            end
        end
    end

    % -- Total readout time (topup) --
    if strcmp(CFG.PREPROC_MODE, 'topup')
        bidsReadout = [];
        if isfield(bj, 'EstimatedTotalReadoutTime')
            bidsReadout = bj.EstimatedTotalReadoutTime;
        elseif isfield(bj, 'TotalReadoutTime')
            bidsReadout = bj.TotalReadoutTime;
        end
        if ~isempty(bidsReadout)
            if abs(CFG.TOPUP_READOUT_SEC - bidsReadout) > 0.001
                fprintf('  WARNING: Readout time mismatch — config TOPUP_READOUT_SEC: %.6f s, BIDS JSON: %.6f s\n', ...
                    CFG.TOPUP_READOUT_SEC, bidsReadout);
                nWarnings = nWarnings + 1;
            else
                fprintf('  Readout time (topup): %.6f s (matches)\n', bidsReadout);
            end
            % Sanity: warn if value looks like milliseconds
            if CFG.TOPUP_READOUT_SEC > 1
                fprintf('  WARNING: TOPUP_READOUT_SEC=%.3f looks like milliseconds — should be in seconds\n', ...
                    CFG.TOPUP_READOUT_SEC);
                nWarnings = nWarnings + 1;
            end
        end
    end

    % -- Total readout time (realign_unwarp — in ms) --
    if strcmp(CFG.PREPROC_MODE, 'realign_unwarp')
        bidsReadout = [];
        if isfield(bj, 'TotalReadoutTime')
            bidsReadout = bj.TotalReadoutTime * 1000;  % BIDS stores seconds
        elseif isfield(bj, 'EstimatedTotalReadoutTime')
            bidsReadout = bj.EstimatedTotalReadoutTime * 1000;
        end
        if ~isempty(bidsReadout)
            if abs(CFG.TOTAL_READOUT - bidsReadout) > 1  % tolerance 1 ms
                fprintf('  WARNING: Readout time mismatch — config TOTAL_READOUT_MS: %.3f ms, BIDS JSON: %.3f ms\n', ...
                    CFG.TOTAL_READOUT, bidsReadout);
                nWarnings = nWarnings + 1;
            else
                fprintf('  Readout time (VDM): %.3f ms (matches)\n', bidsReadout);
            end
        end
    end
else
    fprintf('  No BOLD JSON sidecar found — skipping BOLD parameter validation.\n');
end

% --- The reverse-PE EPI topup was given ---
if strcmp(CFG.PREPROC_MODE, 'topup') && ~isempty(reverseFile)
    revJson = regexprep(reverseFile, '\.nii(\.gz)?$', '.json');

    if isfile(revJson)
        rj = read_json_file(revJson);
        fprintf('  Reverse-PE JSON: %s\n', revJson);

        revPE = '';
        revPE_hasSign = false;
        if isfield(rj, 'PhaseEncodingDirection') && ~isempty(rj.PhaseEncodingDirection)
            revPE = rj.PhaseEncodingDirection;
            revPE_hasSign = true;
        elseif isfield(rj, 'PhaseEncodingAxis') && ~isempty(rj.PhaseEncodingAxis)
            revPE = rj.PhaseEncodingAxis;
            revPE_hasSign = false;
        end
        if ~isempty(revPE)
            fslRevPE = bids_pe_to_fsl(revPE);
            configRevAxis = regexprep(CFG.TOPUP_PE_DIR_REVERSE, '-$', '');
            bidsRevAxis = regexprep(fslRevPE, '-$', '');

            if revPE_hasSign
                if ~isempty(fslRevPE) && ~strcmp(fslRevPE, CFG.TOPUP_PE_DIR_REVERSE)
                    fprintf('  WARNING: Reverse PE direction mismatch — config: %s, BIDS JSON: %s (=%s in FSL)\n', ...
                        CFG.TOPUP_PE_DIR_REVERSE, revPE, fslRevPE);
                    nWarnings = nWarnings + 1;
                else
                    fprintf('  PE direction (reverse): %s = %s in FSL (matches)\n', revPE, fslRevPE);
                end
            else
                if ~strcmp(bidsRevAxis, configRevAxis)
                    fprintf('  WARNING: Reverse PE axis mismatch — config axis: %s, BIDS JSON: %s (=%s in FSL)\n', ...
                        configRevAxis, revPE, bidsRevAxis);
                    nWarnings = nWarnings + 1;
                else
                    fprintf('  PE axis (reverse): %s (matches). PhaseEncodingAxis only — cannot verify +/- sign.\n', bidsRevAxis);
                end
            end
        end

        % Forward and reverse must be opposite, or topup has nothing to work with
        if bidsPE_hasSign && revPE_hasSign && ~isempty(bidsPE) && ~isempty(revPE)
            fwdVec = pe_dir_to_vector(bids_pe_to_fsl(bidsPE));
            revVec = pe_dir_to_vector(bids_pe_to_fsl(revPE));
            if ~all(fwdVec + revVec == 0)
                fprintf('  WARNING: Forward PE (%s) and reverse PE (%s) are NOT opposite!\n', bidsPE, revPE);
                fprintf('           Forward vector: [%d %d %d], Reverse vector: [%d %d %d]\n', ...
                    fwdVec, revVec);
                nWarnings = nWarnings + 1;
            else
                fprintf('  Forward/reverse PE vectors are opposite (OK)\n');
            end
        elseif ~isempty(bidsPE) && ~isempty(revPE)
            fwdAxis = regexprep(bids_pe_to_fsl(bidsPE), '-$', '');
            revAxis = regexprep(bids_pe_to_fsl(revPE), '-$', '');
            if strcmp(fwdAxis, revAxis)
                fprintf('  Forward and reverse PE share axis %s (OK). Only PhaseEncodingAxis available — verify opposite signs manually.\n', fwdAxis);
            else
                fprintf('  WARNING: Forward PE axis (%s) and reverse PE axis (%s) are on DIFFERENT axes!\n', fwdAxis, revAxis);
                nWarnings = nWarnings + 1;
            end
        end

        % The two EPIs must also come from the same visit to the scanner, or
        % topup models the head moving as a distortion field.
        nWarnings = nWarnings + warn_if_unmatched(boldFiles, reverseFile, 'reverse-PE EPI', CFG);
    else
        fprintf('  No reverse-PE JSON sidecar found — skipping reverse-PE validation.\n');
    end
end

% --- The GRE phasediff ---
if strcmp(CFG.PREPROC_MODE, 'realign_unwarp') && ~isempty(phasediffFile)
    phasediffJson = regexprep(phasediffFile, '\.nii(\.gz)?$', '.json');
    if isfile(phasediffJson)
        pj = read_json_file(phasediffJson);
        fprintf('  Phasediff JSON: %s\n', phasediffJson);

        if isfield(pj, 'EchoTime1') && isfield(pj, 'EchoTime2')
            bidsTE1 = pj.EchoTime1 * 1000;  % to ms
            bidsTE2 = pj.EchoTime2 * 1000;
            if abs(CFG.TE_SHORT - bidsTE1) > 0.01
                fprintf('  WARNING: TE1 mismatch — config TE_SHORT_MS: %.3f ms, BIDS JSON EchoTime1: %.3f ms\n', ...
                    CFG.TE_SHORT, bidsTE1);
                nWarnings = nWarnings + 1;
            else
                fprintf('  TE1: %.3f ms (matches)\n', bidsTE1);
            end
            if abs(CFG.TE_LONG - bidsTE2) > 0.01
                fprintf('  WARNING: TE2 mismatch — config TE_LONG_MS: %.3f ms, BIDS JSON EchoTime2: %.3f ms\n', ...
                    CFG.TE_LONG, bidsTE2);
                nWarnings = nWarnings + 1;
            else
                fprintf('  TE2: %.3f ms (matches)\n', bidsTE2);
            end
        else
            fprintf('  Phasediff JSON has no EchoTime1/EchoTime2 — skipping TE validation.\n');
        end
    else
        fprintf('  No phasediff JSON sidecar — skipping TE validation.\n');
    end
    nWarnings = nWarnings + warn_if_unmatched(boldFiles, phasediffFile, 'phasediff', CFG);
end

% --- The precalculated fieldmap ---
if strcmp(CFG.PREPROC_MODE, 'precalc_fieldmap') && ~isempty(fieldmapFile)
    fmJson = regexprep(fieldmapFile, '\.nii(\.gz)?$', '.json');
    if ~isfile(fmJson)
        fprintf('  No fieldmap JSON next to %s — cannot cross-check units.\n', fieldmapFile);
        fprintf('  Using FIELDMAP_UNITS from config: %s\n', CFG.FIELDMAP_UNITS);
    else
        fj = read_json_file(fmJson);
        fprintf('  Fieldmap JSON: %s\n', fmJson);
        if isfield(fj, 'Units') && ~isempty(fj.Units)
            jsonUnits = lower(strtrim(fj.Units));
            cfgUnits  = CFG.FIELDMAP_UNITS;
            % 'rad/s' and 'rads' mean the same thing here
            normalise = @(u) regexprep(lower(u), '^(rad/s|rads|rad_per_s)$', 'rad/s');
            if ~strcmp(normalise(jsonUnits), normalise(cfgUnits))
                fprintf(['  WARNING: fieldmap units mismatch — config FIELDMAP_UNITS: %s, ' ...
                         'JSON Units: %s\n'], cfgUnits, jsonUnits);
                fprintf('           A wrong unit scales the whole distortion correction by 2*pi.\n');
                nWarnings = nWarnings + 1;
            else
                fprintf('  Fieldmap units: %s (matches config)\n', jsonUnits);
            end
        else
            fprintf('  Fieldmap JSON has no "Units" field — using config: %s\n', CFG.FIELDMAP_UNITS);
        end
        if isfield(fj, 'GeneratedBy') && isstruct(fj.GeneratedBy) && ...
           isfield(fj.GeneratedBy, 'MagnitudeSource')
            fprintf('  Magnitude source: %s\n', fj.GeneratedBy.MagnitudeSource);
        end
    end

    % A GRE-derived fieldmap is not EPI-based; getting this backwards inverts
    % how SPM matches the VDM to the EPI.
    if CFG.EPI_BASED_FIELDMAP ~= 0
        fprintf(['  WARNING: EPI_BASED_FIELDMAP=%d with PREPROC_MODE=precalc_fieldmap.\n' ...
                 '           A fieldmap made from a GRE phasediff (make_fieldmaps.sh)\n' ...
                 '           is NOT EPI-based — set EPI_BASED_FIELDMAP=0.\n'], ...
                CFG.EPI_BASED_FIELDMAP);
        nWarnings = nWarnings + 1;
    end
end

% --- Summary ---
if nWarnings == 0
    fprintf('  All BIDS JSON checks passed.\n');
else
    fprintf('\n  *** %d BIDS JSON warning(s) detected. ***\n', nWarnings);
    fprintf('  Config values will still be used. Check warnings above.\n');
    fprintf('  If the config values are intentional overrides, you can ignore these warnings.\n');
end
fprintf('\n');
end

function n = warn_if_unmatched(boldFiles, otherFile, label, CFG)
% WARN_IF_UNMATCHED  Say so when a chosen input was acquired far from the BOLD.
%
% The run selection already tries to match by acquisition time, but an explicit
% choice in run_selection.tsv overrides that, and a missing sidecar hides it.
% This is the last place to notice that a distortion correction is about to be
% estimated between scans from two different head positions.

n = 0;
if isempty(otherFile), return; end

otherSec = bids_acq_time(otherFile);
if isnan(otherSec), return; end

boldSecs = [];
for i = 1:numel(boldFiles)
    if isempty(boldFiles{i}), continue; end
    s = bids_acq_time(boldFiles{i});
    if ~isnan(s), boldSecs(end+1) = s; end %#ok<AGROW>
end
if isempty(boldSecs), return; end

gap = min(abs(boldSecs - otherSec));
if gap > CFG.RUN_MATCH_GAP_SEC
    fprintf(['  WARNING: the %s was acquired %.0f min from the nearest BOLD run.\n' ...
             '           They are probably from different visits to the scanner —\n' ...
             '           check run_selection.tsv.\n'], label, gap/60);
    n = 1;
end
end

function fslDir = bids_pe_to_fsl(bidsPE)
% BIDS_PE_TO_FSL  Convert BIDS PhaseEncodingDirection to FSL convention.
%   'i' -> 'x',  'i-' -> 'x-'
%   'j' -> 'y',  'j-' -> 'y-'
%   'k' -> 'z',  'k-' -> 'z-'
%   Already FSL format (x/y/z) passes through unchanged.
%   Returns '' if unknown.

switch bidsPE
    case 'i',  fslDir = 'x';
    case 'i-', fslDir = 'x-';
    case 'j',  fslDir = 'y';
    case 'j-', fslDir = 'y-';
    case 'k',  fslDir = 'z';
    case 'k-', fslDir = 'z-';
    case {'x','x-','y','y-','z','z-'}
        fslDir = bidsPE;  % Already FSL format
    otherwise
        fslDir = '';
        fprintf('  (Unknown BIDS PE direction: %s)\n', bidsPE);
end
end

% ===================== PROVENANCE LOG =====================

function save_provenance_log(DER, SUB, SES, CFG, cfgFile, taskBoldFiles, boldPicks, stagesRun)
% SAVE_PROVENANCE_LOG  Write a JSON log of all preprocessing settings used.
%
% Saved as preproc_provenance.json at the top of the session's derivatives.
% Contains the preprocessing mode, the stages this run executed, which run of
% each scan it chose, all key parameters, timestamps, software versions, and
% input files.
%
% The JSON always describes the most recent run. Because stages can be run
% separately, a running history of (timestamp, mode, stages) is also appended
% to preproc_stages.log next to it.

funcDir = DER.func;
logFile = fullfile(DER.sesDir, 'preproc_provenance.json');
fprintf('\nSaving provenance log: %s\n', logFile);

prov = struct();

% --- Metadata ---
prov.created = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
prov.subject = SUB;
prov.session = SES;
prov.config_file = cfgFile;
prov.preproc_mode = CFG.PREPROC_MODE;
prov.stages_run = stagesRun;

% --- Software versions ---
prov.software.spm_dir = CFG.SPM_DIR;
try
    prov.software.spm_version = spm('Version');
catch
    prov.software.spm_version = 'unknown';
end
prov.software.matlab_version = version;
if strcmp(CFG.PREPROC_MODE, 'topup')
    [st, fslVer] = system('cat $FSLDIR/etc/fslversion 2>/dev/null');
    if st == 0
        prov.software.fsl_version = strtrim(fslVer);
    else
        prov.software.fsl_version = 'unknown';
    end
end

% --- Input files, and which run of each was chosen ---
prov.inputs.task_labels = CFG.TASK_LABELS;
prov.inputs.bold_files = taskBoldFiles;
prov.inputs.derivatives_dir = DER.sesDir;
prov.inputs.deriv_reset = CFG.DERIV_RESET;

prov.run_selection.file = CFG.RUN_SELECTION_FILE;
prov.run_selection.strict = CFG.STRICT_RUN_MATCHING;
prov.run_selection.match_gap_min = CFG.RUN_MATCH_GAP_SEC / 60;
for i = 1:numel(boldPicks)
    if isempty(boldPicks(i).path), continue; end
    key = matlab.lang.makeValidName(['task_' boldPicks(i).label]);
    prov.run_selection.chosen.(key) = struct( ...
        'run',    boldPicks(i).runLabel, ...
        'source', boldPicks(i).path, ...
        'reason', boldPicks(i).reason);
end

% --- Acquisition parameters ---
if isfield(CFG, 'TR'), prov.acquisition.TR = CFG.TR; end
prov.acquisition.total_readout_ms = CFG.TOTAL_READOUT;
prov.acquisition.te_short_ms = CFG.TE_SHORT;
prov.acquisition.te_long_ms = CFG.TE_LONG;
prov.acquisition.blip_direction = CFG.BLIPDIR;
prov.acquisition.epi_based_fieldmap = CFG.EPI_BASED_FIELDMAP;

% --- Precalculated fieldmap parameters (if applicable) ---
if strcmp(CFG.PREPROC_MODE, 'precalc_fieldmap')
    prov.fieldmap.units = CFG.FIELDMAP_UNITS;
    prov.fieldmap.pattern = CFG.FIELDMAP_PATTERN;
    prov.fieldmap.magnitude_pattern = CFG.FIELDMAP_MAGNITUDE_PATTERN;
end

% --- Topup parameters (if applicable) ---
if strcmp(CFG.PREPROC_MODE, 'topup')
    prov.topup.pe_dir_bold = CFG.TOPUP_PE_DIR_BOLD;
    prov.topup.pe_dir_reverse = CFG.TOPUP_PE_DIR_REVERSE;
    prov.topup.readout_sec = CFG.TOPUP_READOUT_SEC;
    prov.topup.apply_method = CFG.TOPUP_APPLY_METHOD;
    prov.topup.interp = CFG.TOPUP_INTERP;
    prov.topup.reverse_pe_pattern = CFG.TOPUP_REVERSE_PE_PATTERN;
    prov.topup.config = CFG.TOPUP_CONFIG;
end

% --- Realign parameters ---
prov.realign.quality = CFG.REALIGN_QUALITY;
prov.realign.sep = CFG.REALIGN_SEP;
prov.realign.fwhm = CFG.REALIGN_FWHM;
prov.realign.rtm = CFG.REALIGN_RTM;
prov.realign.einterp = CFG.REALIGN_EINTERP;
prov.realign.ewrap = CFG.REALIGN_EWRAP;

% --- Unwarp parameters ---
prov.unwarp.basfcn = CFG.UNWARP_BASFCN;
prov.unwarp.regorder = CFG.UNWARP_REGORDER;
prov.unwarp.lambda = CFG.UNWARP_LAMBDA;
prov.unwarp.jm = CFG.UNWARP_JM;
prov.unwarp.fot = CFG.UNWARP_FOT;
prov.unwarp.sot = CFG.UNWARP_SOT;
prov.unwarp.uwfwhm = CFG.UNWARP_UWFWHM;
prov.unwarp.rem = CFG.UNWARP_REM;
prov.unwarp.noi = CFG.UNWARP_NOI;
prov.unwarp.expround = CFG.UNWARP_EXPROUND;

% --- Reslice parameters ---
prov.reslice.which = CFG.WHICHUNWARP;
prov.reslice.interp = CFG.RESLICE_INTERP;
prov.reslice.wrap = CFG.RESLICE_WRAP;
prov.reslice.mask = CFG.RESLICE_MASK;
prov.reslice.prefix = CFG.RESLICE_PREFIX;

% --- VDM parameters ---
prov.vdm.wrap = CFG.VDM_WRAP;
prov.vdm.mask = CFG.VDM_MASK;
prov.vdm.maskbrain = CFG.MASKBRAIN;
prov.vdm.uflags_method = CFG.VDM_UFLAGS_METHOD;
prov.vdm.uflags_fwhm = CFG.VDM_UFLAGS_FWHM;
prov.vdm.uflags_pad = CFG.VDM_UFLAGS_PAD;
prov.vdm.uflags_ws = CFG.VDM_UFLAGS_WS;

% --- Coregistration parameters ---
prov.coreg.cost_fun = CFG.COREG_COST_FUN;
prov.coreg.sep = CFG.COREG_SEP;
prov.coreg.fwhm = CFG.COREG_FWHM;
prov.coreg.roi_interp = CFG.COREG_ROI_INTERP;
prov.coreg.do_coreg_rois = CFG.DO_COREG_ROIS;

% --- Smoothing parameters ---
if isfield(CFG, 'SMOOTH_FWHM')
    prov.smoothing.fwhm = CFG.SMOOTH_FWHM;
    prov.smoothing.prefix = CFG.SMOOTH_PREFIX;
    prov.smoothing.applied = CFG.SMOOTH_FWHM > 0;
end

% --- BIDS JSON values (for cross-reference) ---
% The sidecar of the BOLD run this session actually used — not whichever
% sidecar happened to sort first, which said the wrong thing for multi-run
% sessions.
boldJson = '';
for i = 1:numel(taskBoldFiles)
    if isempty(taskBoldFiles{i}), continue; end
    cand = regexprep(taskBoldFiles{i}, '\.nii(\.gz)?$', '.json');
    if isfile(cand)
        boldJson = cand;
        break;
    end
end
if ~isempty(boldJson)
    bj = read_json_file(boldJson);
    prov.bids_json.source_file = boldJson;
    if isfield(bj, 'RepetitionTime')
        prov.bids_json.RepetitionTime = bj.RepetitionTime;
    end
    if isfield(bj, 'EstimatedTotalReadoutTime')
        prov.bids_json.EstimatedTotalReadoutTime = bj.EstimatedTotalReadoutTime;
    end
    if isfield(bj, 'TotalReadoutTime')
        prov.bids_json.TotalReadoutTime = bj.TotalReadoutTime;
    end
    if isfield(bj, 'PhaseEncodingDirection')
        prov.bids_json.PhaseEncodingDirection = bj.PhaseEncodingDirection;
    end
    if isfield(bj, 'PhaseEncodingAxis')
        prov.bids_json.PhaseEncodingAxis = bj.PhaseEncodingAxis;
    end
    if isfield(bj, 'EffectiveEchoSpacing')
        prov.bids_json.EffectiveEchoSpacing = bj.EffectiveEchoSpacing;
    end
    if isfield(bj, 'EstimatedEffectiveEchoSpacing')
        prov.bids_json.EstimatedEffectiveEchoSpacing = bj.EstimatedEffectiveEchoSpacing;
    end
    if isfield(bj, 'EchoTime')
        prov.bids_json.EchoTime = bj.EchoTime;
    end
end

% --- Write JSON ---
jsonText = jsonencode(prov);
% Pretty-print: add newlines after { and , for readability
jsonText = strrep(jsonText, ',"', sprintf(',\n  "'));
jsonText = strrep(jsonText, '{"', sprintf('{\n  "'));
jsonText = strrep(jsonText, '}', sprintf('\n}'));

fid = fopen(logFile, 'w');
if fid > 0
    fprintf(fid, '%s\n', jsonText);
    fclose(fid);
    fprintf('  Provenance log saved.\n');
else
    warning('Could not write provenance log: %s', logFile);
end

% --- Append this run to the stage history ---
stageLog = fullfile(funcDir, 'preproc_stages.log');
fid = fopen(stageLog, 'a');
if fid > 0
    fprintf(fid, '%s\t%s\t%s\t%s\n', prov.created, CFG.PREPROC_MODE, ...
            strjoin(stagesRun, ','), cfgFile);
    fclose(fid);
end
end

% ===================== STAGE HELPERS =====================

function S = resolve_stages(spec)
% RESOLVE_STAGES  Turn a stage specification into a struct of flags.
%
% spec may be:
%   ''  or  'all'                     -> every stage
%   'topup,vdm,realign,coreg,smooth'  -> an explicit comma-separated list
%   {'realign','coreg'}               -> the same as a cell array
%   'fieldmap'                        -> shorthand for topup,vdm
%   'post_fieldmap'                   -> shorthand for realign,coreg,smooth
%
% Returns a struct with one logical field per stage, plus:
%   .list     stage names to run, in canonical order
%   .skipped  stage names that will not run, in canonical order

ALL = {'topup','vdm','realign','coreg','smooth'};

if isempty(spec)
    spec = 'all';
end

if isstring(spec) && ~isscalar(spec)
    spec = cellstr(spec);   % a string array behaves like a cellstr here
end

if ischar(spec) || isstring(spec)
    parts = strsplit(char(spec), ',');
elseif iscell(spec)
    parts = spec;
else
    error('run_spm_preproc:BadStageSpec', ...
          'Stage specification must be a string or cell array of strings.');
end

parts = strtrim(lower(cellfun(@char, parts, 'UniformOutput', false)));
parts = parts(~cellfun(@isempty, parts));

sel = {};
for i = 1:numel(parts)
    p = parts{i};
    switch p
        case 'all'
            sel = [sel, ALL]; %#ok<AGROW>
        case {'fieldmap','fieldmaps','fmap'}
            sel = [sel, {'topup','vdm'}]; %#ok<AGROW>
        case {'post_fieldmap','postfieldmap','post-fieldmap'}
            sel = [sel, {'realign','coreg','smooth'}]; %#ok<AGROW>
        case {'unwarp','realign_unwarp'}
            sel = [sel, {'realign'}]; %#ok<AGROW>
        otherwise
            if ~ismember(p, ALL)
                error('run_spm_preproc:UnknownStage', [ ...
                      'Unknown preprocessing stage: "%s"\n' ...
                      'Valid stages: %s\n' ...
                      'Shorthands:   all, fieldmap (topup,vdm), post_fieldmap (realign,coreg,smooth)'], ...
                      p, strjoin(ALL, ', '));
            end
            sel{end+1} = p; %#ok<AGROW>
    end
end

assert(~isempty(sel), 'run_spm_preproc:NoStages', 'No preprocessing stages selected.');

S = struct();
for i = 1:numel(ALL)
    S.(ALL{i}) = ismember(ALL{i}, sel);
end
S.list    = ALL(ismember(ALL, sel));      % canonical order
S.skipped = ALL(~ismember(ALL, sel));
end

function topupPrefix = resolve_topup_prefix(funcDir, SUB, SES, CFG, useExisting)
% RESOLVE_TOPUP_PREFIX  Where the FSL topup output for this session lives.
%
% The topup stage always writes into func/ as "topup_results". When that stage
% is skipped, TOPUP_EXISTING_PREFIX (if set) points at topup output produced
% outside this pipeline; {SUB} and {SES} in it are substituted. A full filename
% such as ".../mytopup_fieldcoef.nii.gz" is accepted and reduced to its prefix.

topupPrefix = fullfile(funcDir, 'topup_results');

if ~useExisting || isempty(CFG.TOPUP_EXISTING_PREFIX)
    return;
end

p = CFG.TOPUP_EXISTING_PREFIX;
p = strrep(p, '{SUB}', SUB);
p = strrep(p, '{SES}', SES);
p = regexprep(p, '_(fieldcoef|field|corrected)\.nii(\.gz)?$', '');
p = regexprep(p, '_movpar\.txt$', '');
topupPrefix = p;
end

function assert_topup_outputs(topupPrefix, CFG, STAGES)
% ASSERT_TOPUP_OUTPUTS  Check that pre-computed topup output is usable.
%
% Only the files the selected stages actually read are required: the field map
% (<prefix>_field) feeds the vdm stage, the spline coefficients and movement
% parameters (<prefix>_fieldcoef, <prefix>_movpar) feed applytopup.

groups = {};
if STAGES.vdm && strcmp(CFG.TOPUP_APPLY_METHOD, 'vdm')
    groups{end+1} = {[topupPrefix '_field.nii.gz'], [topupPrefix '_field.nii']};
end
if STAGES.realign && strcmp(CFG.TOPUP_APPLY_METHOD, 'applytopup')
    groups{end+1} = {[topupPrefix '_fieldcoef.nii.gz'], [topupPrefix '_fieldcoef.nii']};
    groups{end+1} = {[topupPrefix '_movpar.txt']};
end

if isempty(groups)
    fprintf('  No topup output needed by the selected stages.\n');
    return;
end

missing = {};
for i = 1:numel(groups)
    found = '';
    for k = 1:numel(groups{i})
        if isfile(groups{i}{k})
            found = groups{i}{k};
            break;
        end
    end
    if isempty(found)
        missing{end+1} = strjoin(groups{i}, ' or '); %#ok<AGROW>
    else
        fprintf('  Found: %s\n', found);
    end
end

assert(isempty(missing), [ ...
    'Pre-computed FSL topup output is missing:\n  %s\n\n' ...
    'The "topup" stage was not selected, so these files have to exist already.\n' ...
    'Either add "topup" to the stage list, copy your topup output next to\n' ...
    '  %s\n' ...
    'or set TOPUP_EXISTING_PREFIX in pipeline_config.cfg to point at it.'], ...
    strjoin(missing, sprintf('\n  ')), topupPrefix);
end

function write_acqparams(acqFile, CFG)
% WRITE_ACQPARAMS  Write the two-line FSL acqparams file from the config.
% Line 1 = forward (BOLD) phase-encode direction, line 2 = reverse.

fid = fopen(acqFile, 'w');
assert(fid > 0, 'Cannot create acqparams file: %s', acqFile);

fwdVec = pe_dir_to_vector(CFG.TOPUP_PE_DIR_BOLD);
revVec = pe_dir_to_vector(CFG.TOPUP_PE_DIR_REVERSE);
readout = CFG.TOPUP_READOUT_SEC;

fprintf(fid, '%d %d %d %.6f\n', fwdVec(1), fwdVec(2), fwdVec(3), readout);
fprintf(fid, '%d %d %d %.6f\n', revVec(1), revVec(2), revVec(3), readout);
fclose(fid);
end

function acqFile = resolve_acqparams(topupPrefix, funcDir, CFG)
% RESOLVE_ACQPARAMS  Locate the acqparams file applytopup needs.
%
% Prefers the one sitting next to the topup output, since that is what topup
% itself was run with. If topup output produced elsewhere came without one, a
% matching file is written into func/ from TOPUP_PE_DIR_* and TOPUP_READOUT_SEC.

acqFile = fullfile(fileparts(topupPrefix), 'topup_acqparams.txt');
if isfile(acqFile)
    return;
end

acqFile = fullfile(funcDir, 'topup_acqparams.txt');
if ~isfile(acqFile)
    fprintf('  No acqparams file next to the topup output — writing %s from config.\n', acqFile);
    write_acqparams(acqFile, CFG);
end
end

function vdmPath = existing_task_vdm(funcDir, taskLabel, taskIdx)
% EXISTING_TASK_VDM  Path to a voxel displacement map already on disk, or ''.
%
% Matches what calc_vdm_for_run_to_func / convert_topup_to_vdm write
% (vdm_task-<label>.nii). The index form vdm_task-<N>.nii is still accepted:
% it is what earlier versions wrote, and what the documented "drop in your own
% VDMs and run realign,coreg,smooth" workflow tells people to create.

vdmPath = '';
candidates = { fullfile(funcDir, sprintf('vdm_task-%s.nii', taskLabel)), ...
               fullfile(funcDir, sprintf('vdm_task-%s.img', taskLabel)), ...
               fullfile(funcDir, sprintf('vdm_task-%d.nii', taskIdx)), ...
               fullfile(funcDir, sprintf('vdm_task-%d.img', taskIdx)) };
for k = 1:numel(candidates)
    if isfile(candidates{k})
        vdmPath = candidates{k};
        return;
    end
end
end

function sessionMean = resolve_session_mean(funcDir, CFG, rebuild)
% RESOLVE_SESSION_MEAN  Find (or build) the mean EPI used as coreg reference.
%
%   rebuild = true   always recompute from the per-task means; used right after
%                    realignment, when those means have just been rewritten
%   rebuild = false  reuse an existing session mean if there is one; used when
%                    the coreg stage runs on its own

sessionMeanFile = fullfile(funcDir, CFG.SESSION_MEAN_NAME);

if ~rebuild && isfile(sessionMeanFile)
    sessionMean = sessionMeanFile;
    fprintf('\nUsing existing session mean: %s\n', sessionMean);
    return;
end

% --------- Collect the mean image written for each task ----------
% SPM Realign&Unwarp writes "meanu*" while Realign-only writes "mean*".
% Try task-specific match first, then broad fallback.
meanu_imgs = cell(0,1);
for ii = 1:numel(CFG.TASK_LABELS)
    tLabel = CFG.TASK_LABELS{ii};
    escLabel = regexptranslate('escape', tLabel);
    % The boundary after the label matters: without it task-run1 also matches
    % the mean image of task-run10.
    m = newest_match(funcDir, ['^meanu.*_task-' escLabel '(_|\.).*\.nii$']);
    if isempty(m)
        m = newest_match(funcDir, ['^mean.*_task-' escLabel '(_|\.).*\.nii$']);
    end
    if isempty(m)
        % Broad fallback: any mean* matching this task
        m = newest_match(funcDir, ['^mean.*' escLabel '.*\.nii$']);
    end
    if ~isempty(m)
        meanu_imgs{end+1,1} = m; %#ok<AGROW>
        fprintf('Found mean image (task-%s): %s\n', tLabel, m);
    else
        warning('Could not find mean image for task-%s in %s', tLabel, funcDir);
    end
end

% SPM Realign-only produces a single mean from the first session, not one
% per task.  Deduplicate so the same file isn't counted twice.
meanu_imgs = unique(meanu_imgs);

if CFG.DO_SESSION_MEAN && numel(meanu_imgs) >= 2
    % Multiple per-task means: average them into a session mean
    sessionMean = sessionMeanFile;
    make_session_mean(meanu_imgs, sessionMean);
elseif ~isempty(meanu_imgs)
    % Only one mean available (e.g. Realign-only produces a single mean):
    % use it directly as the session mean / coreg reference.
    sessionMean = meanu_imgs{1};
    fprintf('\nUsing single mean as coreg ref: %s\n', sessionMean);
else
    % No per-task mean at all — try broad fallback
    sessionMean = newest_match(funcDir, '^mean.*\.nii$');
    if isempty(sessionMean)
        warning('No mean image found in %s. Run the "realign" stage first.', funcDir);
    else
        fprintf('\nUsing fallback mean as coreg ref: %s\n', sessionMean);
    end
end
end
