function run_spm_preproc(BIDS_ROOT, SUB, SES, cfgFile, stages)
% run_spm_preproc(BIDS_ROOT, SUB, SES, cfgFile, stages)
%
% Example:
%   run_spm_preproc('/sc-projects/.../b2_bids','sub-01s06','ses-01')
%   run_spm_preproc('/sc-projects/.../b2_bids','sub-01s06','ses-01', '/path/to/config.sh')
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
% from. See DERIV_RESET in config.sh.
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
% STAGES (5th argument, or PREPROC_STAGES in config.sh)
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
%            Set COREG_VIA_TASK to route this through a whole-brain EPI when
%            the T2w and the task EPI are both partial slabs; the masks are
%            still resampled only once. Fails rather than writing a mask with
%            fewer than COREG_MIN_ROI_VOXELS voxels.
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
% Settings come from config.sh (paths, modes) and spm_options.m (SPM parameters).

% --------- Handle optional arguments ----------
if nargin < 4, cfgFile = ''; end
if nargin < 5, stages  = ''; end

if isempty(cfgFile)
    % Default to config.sh next to this script
    scriptDir = fileparts(mfilename('fullpath'));
    cfgFile = fullfile(scriptDir, 'config.sh');
end

assert(isfile(cfgFile), 'config.sh not found at: %s\nProvide cfgFile as 4th argument to run_spm_preproc.', cfgFile);

% spm_options.m is the sibling of config.sh: the two together ARE the
% settings. Putting that directory first on the path means a snapshotted
% run reads the snapshot's settings, not whatever the live checkout has.
addpath(fileparts(cfgFile), '-begin');

% --------- Read configuration ----------
% Two files, each read by the language that uses it:
%   config.sh       paths, modes and filename patterns, as plain strings
%   spm_options.m   everything that goes into an SPM batch, as MATLAB values
%
% There is no mapping table. Anything SPM accepts can be added to
% spm_options.m and it reaches the batch without a line of code here.
cfgRaw = read_config(cfgFile);
OPT    = spm_options();

CFG     = struct();
CFG.spm = OPT;          % the batch builders hand these to SPM wholesale

% ---- from config.sh ----
CFG.SPM_DIR           = cfgRaw.SPM_DIR;
CFG.PREPROC_MODE      = cfgRaw.PREPROC_MODE;
CFG.TASK_LABELS       = strtrim(strsplit(cfgRaw.TASKS, ','));
CFG.DO_COREG_ROIS     = str2double(cfg_or(cfgRaw, 'DO_COREG_ROIS', '1'));
% Empty = coregister the T2w straight to the session mean. A task label
% (usually 'rest') routes it through that whole-brain EPI instead.
CFG.COREG_VIA_TASK    = strtrim(cfg_or(cfgRaw, 'COREG_VIA_TASK', ''));
CFG.COREG_MIN_ROI_VOXELS = str2double(cfg_or(cfgRaw, 'COREG_MIN_ROI_VOXELS', '5'));
CFG.ROI_PATTERN_LEFT  = cfgRaw.ROI_PATTERN_LEFT;
CFG.ROI_PATTERN_RIGHT = cfgRaw.ROI_PATTERN_RIGHT;

% Which run to take when a session holds several of the same scan.
CFG.RUN_SELECTION_FILE  = cfg_or(cfgRaw, 'RUN_SELECTION_FILE', '');
CFG.RUN_MATCH_GAP_SEC   = str2double(cfg_or(cfgRaw, 'RUN_MATCH_GAP_MIN', '20')) * 60;
CFG.STRICT_RUN_MATCHING = cfg_flag(cfgRaw, 'STRICT_RUN_MATCHING', true);

% Clearing the derivatives before a full run: auto | always | never
CFG.DERIV_RESET = lower(cfg_or(cfgRaw, 'DERIV_RESET', 'auto'));

% Finding the fieldmap inputs. The search appends the extension, so strip
% one if the pattern was written with it (the ROI patterns do carry theirs).
stripExt = @(x) regexprep(x, '\.nii(\.gz)?$', '');
CFG.FIELDMAP_PATTERN           = stripExt(cfg_or(cfgRaw, 'FIELDMAP_PATTERN', '_fieldmap'));
CFG.FIELDMAP_MAGNITUDE_PATTERN = stripExt(cfg_or(cfgRaw, 'FIELDMAP_MAGNITUDE_PATTERN', '_magnitude'));
CFG.FIELDMAP_UNITS             = lower(cfg_or(cfgRaw, 'FIELDMAP_UNITS', 'rad/s'));

CFG.TOPUP_APPLY_METHOD       = cfg_or(cfgRaw, 'TOPUP_APPLY_METHOD', 'vdm');
CFG.TOPUP_REVERSE_PE_DIR     = cfg_or(cfgRaw, 'TOPUP_REVERSE_PE_DIR', 'auto');
CFG.TOPUP_REVERSE_PE_PATTERN = cfg_or(cfgRaw, 'TOPUP_REVERSE_PE_PATTERN', '');
CFG.TOPUP_CONFIG             = cfg_or(cfgRaw, 'TOPUP_CONFIG', '');
CFG.TOPUP_EXISTING_PREFIX    = cfg_or(cfgRaw, 'TOPUP_EXISTING_PREFIX', '');

% Output prefixes. config.sh owns these because the shell scripts need them
% too. FUNC_PREFIX is not stored anywhere: it is whatever these two produce.
CFG.RESLICE_PREFIX   = cfg_or(cfgRaw, 'RESLICE_PREFIX',   'u');
CFG.SMOOTH_PREFIX    = cfg_or(cfgRaw, 'SMOOTH_PREFIX',    's');
CFG.ROI_PREFIX       = cfg_or(cfgRaw, 'ROI_PREFIX',       'r');
CFG.SLICETIME_PREFIX = cfg_or(cfgRaw, 'SLICETIME_PREFIX', 'a');

% Slice timing is optional. When it is on, every later stage reads the
% corrected files, so the flag is checked in more than one place.
CFG.DO_SLICE_TIMING = cfg_flag(cfgRaw, 'DO_SLICE_TIMING', false);
CFG.SMOOTH_FWHM    = str2double(cfg_or(cfgRaw, 'SMOOTH_FWHM', '0'));

% ---- from spm_options.m ----
% Batch options are read straight off CFG.spm where they are used. These
% few are also needed outside a batch, so they get a short name here.
CFG.TR                   = OPT.TR;
CFG.DO_SESSION_MEAN      = OPT.session_mean.enabled;
CFG.SESSION_MEAN_NAME    = OPT.session_mean.name;
CFG.TOPUP_PE_DIR_BOLD    = OPT.topup.pe_dir_bold;
CFG.TOPUP_PE_DIR_REVERSE = OPT.topup.pe_dir_reverse;
CFG.TOPUP_READOUT_SEC    = OPT.topup.readout_sec;
CFG.TOPUP_INTERP         = OPT.topup.interp;

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

needVols = STAGES.vdm || STAGES.slicetime || STAGES.realign;

% Smoothing on its own reads only what the earlier stages left in the
% derivatives, so there is nothing to choose and nothing to copy.
%
% Two different questions, and they have different answers for coreg:
%   needSelection  do we have to CHOOSE runs? coreg does - the anchor decides
%                  which T2w, ROI masks and bridge EPI go with these runs.
%   needBoldFiles  do we have to STAGE the BOLD itself? coreg does not: it
%                  reads the session mean that realignment already left in the
%                  derivatives. Staging three decompressed runs it never opens
%                  is the slowest part of a coreg-only re-run.
needSelection = STAGES.topup || STAGES.vdm || STAGES.slicetime || ...
                STAGES.realign || STAGES.coreg;
needBoldFiles = STAGES.topup || STAGES.vdm || STAGES.slicetime || STAGES.realign;

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

    if needBoldFiles
        for ii = 1:nTasks
            if isempty(boldPicks(ii).path), continue; end
            taskBoldFiles{ii} = bids_stage(boldPicks(ii).path, funcDir);
        end
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
bridge = ''; bridgeSrc = '';

if STAGES.coreg && CFG.DO_COREG_ROIS
    t2wSrc = select_t2w(RAW.anat, SUB, SES, CFG, runSel, anchor);
    t2w = bids_stage(t2wSrc, anatDir);
    fprintf('T2w: %s\n', t2w);

    % Optional whole-brain EPI to coregister through. Staged like everything
    % else, because the second estimate rewrites its header.
    if ~isempty(CFG.COREG_VIA_TASK)
        bridgeSrc = select_coreg_bridge(RAW, SUB, SES, CFG, runSel, anchor);
        bridge = bids_stage(bridgeSrc, funcDir);
        fprintf('Coreg bridge: %s\n', bridge);
    end

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
               'reverse', reverseSrc, 'roi_left', roiLsrc, 'roi_right', roiRsrc, ...
               'coreg_bridge', bridgeSrc), ...
        anchor);
end

% --------- Validate config against the BIDS JSON sidecars we chose ----------
% Reading the sidecars of the selected files (not "any file in the folder")
% means the check describes the run actually being preprocessed.
if STAGES.topup || STAGES.vdm || STAGES.realign
    check_bids_params(taskBoldFiles, reversePE_epi, phasemap, fieldmapFile, CFG);
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
                    taskVDM{ii} = make_vdm( ...
                        struct('kind','presubphasemag','phase',phasemap,'magnitude',magnitude1), ...
                        epiRef, funcDir, fmapDir, CFG, tLabel);
                case 'precalc_fieldmap'
                    fprintf('Calculating VDM from precalculated fieldmap for task-%s...\n', tLabel);
                    taskVDM{ii} = make_vdm( ...
                        struct('kind','precalcfieldmap','fieldmap',fieldmapHz,'magnitude',fieldmapMag), ...
                        epiRef, funcDir, fmapDir, CFG, tLabel);
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
%  STAGE slicetime — shift every slice to the same point in the TR
% =====================================================================
%  Optional, and off unless DO_SLICE_TIMING=1 in config.sh. It sits here
%  rather than before the fieldmap stages because the VDM is matched to an
%  untouched EPI reference; slice timing resamples in time only, so the
%  geometry the VDM was built against is unchanged.
%
%  Whether or not this run performs the stage, anything downstream has to
%  read the corrected files if an earlier run wrote them - otherwise
%  resuming with --stages realign,coreg,smooth would silently go back to
%  the uncorrected data.

if STAGES.slicetime
    if ~CFG.DO_SLICE_TIMING
        fprintf('\n=== STAGE slicetime: skipped (DO_SLICE_TIMING=0) ===\n');
    else
        fprintf('\n=== STAGE slicetime: correcting to %.2f x TR ===\n', ...
                CFG.spm.slicetime.ref_fraction);
        corrected = slice_time_correct(taskBoldFiles, CFG);
        for ii = 1:nTasks
            if isempty(corrected{ii}), continue; end
            taskBoldFiles{ii} = corrected{ii};
            taskVols{ii}      = expand_4d(corrected{ii});
        end
    end
elseif CFG.DO_SLICE_TIMING && (STAGES.realign || STAGES.smooth)
    % Resuming after a separate slicetime run: pick up what it wrote.
    for ii = 1:nTasks
        if isempty(taskBoldFiles{ii}), continue; end
        [d, n, e] = fileparts(taskBoldFiles{ii});
        cand = fullfile(d, [CFG.SLICETIME_PREFIX n e]);
        if isfile(cand)
            fprintf('Using slice-time corrected BOLD: %s\n', [CFG.SLICETIME_PREFIX n e]);
            taskBoldFiles{ii} = cand;
            if needVols, taskVols{ii} = expand_4d(cand); end
        elseif STAGES.realign
            error('run_spm_preproc:NoSliceTimedBold', [ ...
                  'DO_SLICE_TIMING=1 but no slice-time corrected BOLD for task-%s:\n' ...
                  '  expected %s\n' ...
                  'Run the "slicetime" stage first, or set DO_SLICE_TIMING=0.'], ...
                  CFG.TASK_LABELS{ii}, cand);
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
            coreg_reslice_rois(sessionMean, t2w, roiList, CFG, bridge);
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

function doReset = deriv_reset_wanted(CFG, STAGES)
% DERIV_RESET_WANTED  Whether to clear this session's derivatives first.
%
% NOTE: this rule also exists in bash, in deriv_reset_wanted() inside
% bids_common.sh. Both are needed and neither is redundant: the bash one
% runs once before an array is submitted (so sessions dropped from the
% subject list do not survive), this one runs per session at job start (so
% calling 01_preprocess.sbatch directly still clears). Change one, change
% the other.
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


function save_provenance_log(DER, SUB, SES, CFG, cfgFile, taskBoldFiles, boldPicks, stagesRun)
% SAVE_PROVENANCE_LOG  Record what this run did, and what it was told to do.
%
% Writes three things to the top of the session's derivatives:
%
%   config.sh          \  verbatim copies of the settings this run used, so
%   spm_options.m      /  the record cannot drift from the code
%   preproc_provenance.json   what ran, on which inputs, with which software
%
% Copying the settings files is deliberate. An earlier version listed every
% parameter by hand, which meant adding a setting took three edits and the
% log silently went stale whenever someone forgot one.
%
% The JSON describes the most recent run. Because stages can be run
% separately, a running history is also appended to preproc_stages.log.

logFile = fullfile(DER.sesDir, 'preproc_provenance.json');
fprintf('\nSaving provenance: %s\n', logFile);

% ---- Copy the settings files next to the output ----
settingsDir = fileparts(cfgFile);
for f = {'config.sh', 'spm_options.m'}
    src = fullfile(settingsDir, f{1});
    if isfile(src)
        copyfile(src, fullfile(DER.sesDir, f{1}));
    else
        warning('Could not copy %s for provenance (not found at %s)', f{1}, src);
    end
end

% ---- What ran ----
prov = struct();
prov.created      = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
prov.subject      = SUB;
prov.session      = SES;
prov.preproc_mode = CFG.PREPROC_MODE;
prov.stages_run   = stagesRun;
prov.settings     = {'config.sh', 'spm_options.m'};   % copied alongside

% ---- Software ----
prov.software.spm_dir        = CFG.SPM_DIR;
prov.software.matlab_version = version;
try
    prov.software.spm_version = spm('Version');
catch
    prov.software.spm_version = 'unknown';
end
if strcmp(CFG.PREPROC_MODE, 'topup')
    [st, fslVer] = system('cat $FSLDIR/etc/fslversion 2>/dev/null');
    prov.software.fsl_version = tern(st == 0, strtrim(fslVer), 'unknown');
end

% ---- Inputs, and which run of each scan was chosen ----
prov.inputs.task_labels     = CFG.TASK_LABELS;
prov.inputs.bold_files      = taskBoldFiles;
prov.inputs.derivatives_dir = DER.sesDir;
prov.inputs.deriv_reset     = CFG.DERIV_RESET;

prov.run_selection.file          = CFG.RUN_SELECTION_FILE;
prov.run_selection.strict        = CFG.STRICT_RUN_MATCHING;
prov.run_selection.match_gap_min = CFG.RUN_MATCH_GAP_SEC / 60;
for i = 1:numel(boldPicks)
    if isempty(boldPicks(i).path), continue; end
    key = matlab.lang.makeValidName(['task_' boldPicks(i).label]);
    prov.run_selection.chosen.(key) = struct( ...
        'run',    boldPicks(i).runLabel, ...
        'source', boldPicks(i).path, ...
        'reason', boldPicks(i).reason);
end

% ---- Write ----
fid = fopen(logFile, 'w');
if fid > 0
    fprintf(fid, '%s\n', jsonencode(prov, 'PrettyPrint', true));
    fclose(fid);
else
    warning('Could not write provenance log: %s', logFile);
end

% Running history, so a resumed run does not erase what came before.
stageLog = fullfile(DER.sesDir, 'preproc_stages.log');
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
%   'topup,vdm,slicetime,realign,coreg,smooth' -> an explicit list
%   {'realign','coreg'}               -> the same as a cell array
%   'fieldmap'                        -> shorthand for topup,vdm
%   'post_fieldmap'                   -> shorthand for slicetime,realign,coreg,smooth
%
% Returns a struct with one logical field per stage, plus:
%   .list     stage names to run, in canonical order
%   .skipped  stage names that will not run, in canonical order

ALL = {'topup','vdm','slicetime','realign','coreg','smooth'};

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
            sel = [sel, {'slicetime','realign','coreg','smooth'}]; %#ok<AGROW>
        case {'stc','slice_timing','slice-timing'}
            sel = [sel, {'slicetime'}]; %#ok<AGROW>
        case {'unwarp','realign_unwarp'}
            sel = [sel, {'realign'}]; %#ok<AGROW>
        otherwise
            if ~ismember(p, ALL)
                error('run_spm_preproc:UnknownStage', [ ...
                      'Unknown preprocessing stage: "%s"\n' ...
                      'Valid stages: %s\n' ...
                      'Shorthands:   all, fieldmap (topup,vdm), post_fieldmap (slicetime,realign,coreg,smooth), stc'], ...
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
    'or set TOPUP_EXISTING_PREFIX in config.sh to point at it.'], ...
    strjoin(missing, sprintf('\n  ')), topupPrefix);
end

function vdmPath = existing_task_vdm(funcDir, taskLabel, taskIdx)
% EXISTING_TASK_VDM  Path to a voxel displacement map already on disk, or ''.
%
% Matches what make_vdm / convert_topup_to_vdm write
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
