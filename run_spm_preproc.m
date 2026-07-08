function run_spm_preproc(BIDS_ROOT, SUB, SES, cfgFile)
% run_spm_preproc(BIDS_ROOT, SUB, SES, cfgFile)
% Example:
%   run_spm_preproc('/sc-projects/.../b2_bids','sub-01s06','ses-01')
%   run_spm_preproc('/sc-projects/.../b2_bids','sub-01s06','ses-01', '/path/to/pipeline_config.cfg')
%
% For one subject/session:
% - Finds fmap magnitude1 + phasediff (1 per session)
% - For each task-1..3:
%     * Calculates run-specific VDM using task's first volume as EPI reference
%       (same phasediff+magnitude1, but matched to each run geometry)
%     * Ensures VDM is moved/renamed into func/ as vdm_task-<n>.nii
% - Runs Realign & Unwarp per session with 3 data blocks
%     (or Realign only, depending on PREPROC_MODE in config)
% - Creates a session mean by averaging meanu(task-1..3)
% - Coregisters T2w->session mean and reslices ROI left/right into EPI space
%
% All settings come from pipeline_config.cfg via read_pipeline_config.m

% --------- Handle optional cfgFile argument ----------
if nargin < 4 || isempty(cfgFile)
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

% Run selection file (for multi-run T2w/BOLD)
if isfield(cfgRaw, 'RUN_SELECTION_FILE') && ~isempty(cfgRaw.RUN_SELECTION_FILE)
    CFG.RUN_SELECTION_FILE = cfgRaw.RUN_SELECTION_FILE;
else
    CFG.RUN_SELECTION_FILE = '';
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

% --------- Safety / init ----------
assert(isfolder(BIDS_ROOT), 'BIDS_ROOT not found: %s', BIDS_ROOT);
assert(startsWith(SUB,'sub-') && startsWith(SES,'ses-'), 'SUB/SES must look like sub-01s06 / ses-01');

addpath(CFG.SPM_DIR);
spm('defaults','FMRI');
spm_jobman('initcfg');

subDir  = fullfile(BIDS_ROOT, SUB, SES);
funcDir = fullfile(subDir, 'func');
anatDir = fullfile(subDir, 'anat');
fmapDir = fullfile(subDir, 'fmap');

assert(isfolder(funcDir), 'Missing func dir: %s', funcDir);
assert(isfolder(anatDir), 'Missing anat dir: %s', anatDir);
if strcmp(CFG.PREPROC_MODE, 'realign_unwarp')
    assert(isfolder(fmapDir), 'Missing fmap dir: %s', fmapDir);
elseif strcmp(CFG.PREPROC_MODE, 'topup')
    % For topup: fmap required unless searching func only
    if strcmp(CFG.TOPUP_REVERSE_PE_DIR, 'fmap')
        assert(isfolder(fmapDir), 'Missing fmap dir: %s (TOPUP_REVERSE_PE_DIR=fmap)', fmapDir);
    end
    % 'auto' and 'func' don't strictly require fmap/
end

fprintf('\n=== %s / %s ===\n', SUB, SES);
fprintf('Preprocessing mode: %s\n', CFG.PREPROC_MODE);

% --------- Load run selection (multi-run overrides) ----------
runSel = load_run_selection(CFG.RUN_SELECTION_FILE, SUB, SES);

% --------- Find fieldmap files ----------
phasemap = '';
magnitude1 = '';
reversePE_epi = '';

if strcmp(CFG.PREPROC_MODE, 'realign_unwarp')
    % GRE fieldmap mode: need phasediff + magnitude1
    phasemap = first_match(fmapDir, [SUB '_' SES '.*phasediff.*\.nii$']);
    if isempty(phasemap)
        phasemap = first_match(fmapDir, '.*phasediff.*\.nii$');
    end
    magnitude1 = first_match(fmapDir, [SUB '_' SES '.*magnitude1.*\.nii$']);
    if isempty(magnitude1)
        magnitude1 = first_match(fmapDir, '.*magnitude1.*\.nii$');
    end

    assert(~isempty(phasemap),   'Could not find phasediff/phasemap in %s', fmapDir);
    assert(~isempty(magnitude1), 'Could not find magnitude1 in %s', fmapDir);

    phasemap = ensure_nii(phasemap);      % decompress .nii.gz for SPM
    magnitude1 = ensure_nii(magnitude1);

    fprintf('Fieldmap phasediff:  %s\n', phasemap);
    fprintf('Fieldmap magnitude1: %s\n', magnitude1);

elseif strcmp(CFG.PREPROC_MODE, 'topup')
    % Topup mode: find the reverse-PE EPI.
    % Search directories based on TOPUP_REVERSE_PE_DIR setting:
    %   'fmap' = fmap/ only, 'func' = func/ only, 'auto' = fmap/ then func/
    reversePE_epi = find_reverse_pe_epi(subDir, SUB, SES, CFG, runSel);
    fprintf('Reverse-PE EPI:  %s\n', reversePE_epi);
    fprintf('Topup apply method: %s\n', CFG.TOPUP_APPLY_METHOD);
end

% --------- Find T2w (handles multi-run) ----------
t2w = find_t2w(anatDir, SUB, SES, runSel);
t2w = ensure_nii(t2w);  % decompress .nii.gz for SPM
fprintf('T2w: %s\n', t2w);

% --------- Find ROIs (flexible matching) ----------
[roiL, roiR, hasRoiL, hasRoiR] = find_rois(anatDir, SUB, SES, CFG);
if hasRoiL, roiL = ensure_nii(roiL); end  % decompress .nii.gz for SPM
if hasRoiR, roiR = ensure_nii(roiR); end
if CFG.DO_COREG_ROIS
    fprintf('ROI left exists:  %d', hasRoiL);
    if hasRoiL, fprintf(' (%s)', roiL); end
    fprintf('\n');
    fprintf('ROI right exists: %d', hasRoiR);
    if hasRoiR, fprintf(' (%s)', roiR); end
    fprintf('\n');
end

% --------- Validate config against BIDS JSON sidecars ----------
validate_bids_params(funcDir, fmapDir, SUB, SES, CFG);

% --------- Prepare per-task scans + per-task VDMs ----------
nTasks = numel(CFG.TASK_LABELS);
taskVols      = cell(nTasks,1);
taskVDM       = cell(nTasks,1);
taskBoldFiles = cell(nTasks,1);
meanu_imgs    = cell(0,1);

for ii = 1:nTasks
    tLabel = CFG.TASK_LABELS{ii};   % e.g. 'run1', '1', 'run2', etc.

    boldFile = find_bold(funcDir, SUB, SES, tLabel, runSel);
    boldFile = ensure_nii(boldFile);  % decompress .nii.gz for SPM

    taskBoldFiles{ii} = boldFile;
    vols = expand_4d(boldFile);     % cellstr of '...nii,1' '...nii,2' ...
    epiRef = vols{1};               % char: first volume as VDM reference

    fprintf('\n--- Task %s ---\n', tLabel);
    fprintf('BOLD: %s\n', boldFile);

    % For realign_unwarp mode: calculate VDM matched to this run's EPI reference
    if strcmp(CFG.PREPROC_MODE, 'realign_unwarp')
        fprintf('Calculating VDM for task-%s...\n', tLabel);
        vdm_out = calc_vdm_for_run_to_func(phasemap, magnitude1, epiRef, funcDir, fmapDir, CFG, ii);
        fprintf('VDM saved: %s\n', vdm_out);
        taskVDM{ii} = vdm_out;
    else
        taskVDM{ii} = '';
    end

    taskVols{ii} = vols;
end

% --------- Run distortion correction + motion correction ----------
fprintf('\nRunning preprocessing job for all tasks...\n');

if strcmp(CFG.PREPROC_MODE, 'realign_unwarp')
    run_realign_unwarp_multi(taskVols, taskVDM, CFG);

elseif strcmp(CFG.PREPROC_MODE, 'topup')
    % --- FSL TOPUP pathway ---
    % Step 1: Run FSL topup to estimate the distortion field (once per session)
    fprintf('\n--- Running FSL topup (session-level) ---\n');
    topupPrefix = fullfile(funcDir, 'topup_results');
    run_fsl_topup(taskBoldFiles{1}, reversePE_epi, topupPrefix, CFG);

    if strcmp(CFG.TOPUP_APPLY_METHOD, 'applytopup')
        % Step 2a: Apply topup correction to each task's 4D BOLD, then realign
        fprintf('\n--- Applying topup correction via applytopup ---\n');
        correctedBoldFiles = cell(nTasks,1);
        for ii = 1:nTasks
            tLabel = CFG.TASK_LABELS{ii};
            correctedBold = fullfile(funcDir, sprintf('%s_%s_task-%s_bold_dc.nii', SUB, SES, tLabel));
            correctedBold = run_fsl_applytopup(taskBoldFiles{ii}, topupPrefix, correctedBold, 1, CFG);
            correctedBoldFiles{ii} = correctedBold;
            fprintf('Corrected BOLD (task-%s): %s\n', tLabel, correctedBold);
        end

        % Re-expand corrected 4D files for SPM realign
        correctedVols = cell(nTasks,1);
        for ii = 1:nTasks
            correctedVols{ii} = expand_4d(correctedBoldFiles{ii});
        end

        % Realign only (distortion already corrected by applytopup)
        run_realign_only_multi(correctedVols, CFG);

    elseif strcmp(CFG.TOPUP_APPLY_METHOD, 'vdm')
        % Step 2b: Convert topup field to SPM VDM, then Realign&Unwarp
        fprintf('\n--- Converting topup fieldmap to SPM VDM ---\n');
        for ii = 1:nTasks
            tLabel = CFG.TASK_LABELS{ii};
            epiRef = taskVols{ii}{1};
            vdm_out = convert_topup_to_vdm(topupPrefix, epiRef, funcDir, CFG, ii);
            taskVDM{ii} = vdm_out;
            fprintf('VDM from topup (task-%s): %s\n', tLabel, vdm_out);
        end
        run_realign_unwarp_multi(taskVols, taskVDM, CFG);
    else
        error('Unknown TOPUP_APPLY_METHOD: %s (must be applytopup or vdm)', CFG.TOPUP_APPLY_METHOD);
    end

else
    % realign_only mode
    run_realign_only_multi(taskVols, CFG);
end

% --------- Collect mean images for each task after processing ----------
% SPM Realign&Unwarp writes "meanu*" while Realign-only writes "mean*".
% Try task-specific match first, then broad fallback.
for ii = 1:nTasks
    tLabel = CFG.TASK_LABELS{ii};
    escLabel = regexptranslate('escape', tLabel);
    % Try meanu* first (Realign&Unwarp), then mean* (Realign-only)
    m = newest_match(funcDir, ['^meanu.*task-' escLabel '.*\.nii$']);
    if isempty(m)
        m = newest_match(funcDir, ['^mean.*task-' escLabel '.*\.nii$']);
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

% --------- Make session mean from per-task means ----------
sessionMean = '';

if CFG.DO_SESSION_MEAN && numel(meanu_imgs) >= 2
    % Multiple per-task means: average them into a session mean
    sessionMean = fullfile(funcDir, CFG.SESSION_MEAN_NAME);
    make_session_mean(meanu_imgs, sessionMean);
elseif ~isempty(meanu_imgs)
    % Only one mean available (e.g. Realign-only produces a single mean):
    % use it directly as the session mean / coreg reference.
    sessionMean = meanu_imgs{1};
    fprintf('\nUsing single mean as coreg ref: %s\n', sessionMean);
else
    % No mean at all — try broad fallback
    sessionMean = newest_match(funcDir, '^mean.*\.nii$');
    if isempty(sessionMean)
        warning('No mean image found in %s. Skipping coreg.', funcDir);
    else
        fprintf('\nUsing fallback mean as coreg ref: %s\n', sessionMean);
    end
end

% --------- Coreg + reslice ROIs into EPI space ----------
if CFG.DO_COREG_ROIS && ~isempty(sessionMean) && (hasRoiL || hasRoiR)
    fprintf('\nCoreg+reslice ROIs into EPI space using session mean...\n');
    roiList = {};
    if hasRoiL, roiList{end+1} = roiL; end %#ok<AGROW>
    if hasRoiR, roiList{end+1} = roiR; end %#ok<AGROW>
    coreg_reslice_rois(sessionMean, t2w, roiList, CFG);
else
    fprintf('\nSkipping ROI coreg (missing session mean or ROIs).\n');
end


% --------- Optional spatial smoothing ----------
if CFG.SMOOTH_FWHM > 0
    fprintf('\nSmoothing preprocessed BOLD files (FWHM = %g mm)...\n', CFG.SMOOTH_FWHM);
    smooth_bold_files(funcDir, CFG);
else
    fprintf('\nSmoothing: skipped (SMOOTH_FWHM = 0).\n');
end
% --------- Save preprocessing provenance log ----------
save_provenance_log(funcDir, SUB, SES, CFG, cfgFile, taskBoldFiles);

fprintf('\nDONE: %s / %s\n\n', SUB, SES);
end

% ============================== HELPERS ==============================

function out = first_match(folder, regexPattern)
d = dir(folder);
out = '';
for i=1:numel(d)
    if d(i).isdir, continue; end
    if ~isempty(regexp(d(i).name, regexPattern, 'once'))
        out = fullfile(folder, d(i).name);
        return;
    end
end
end

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

function matches = all_matches(folder, regexPattern)
% Return all files matching regex in folder (cell array of full paths)
d = dir(folder);
matches = {};
for i = 1:numel(d)
    if d(i).isdir, continue; end
    if ~isempty(regexp(d(i).name, regexPattern, 'once'))
        matches{end+1} = fullfile(folder, d(i).name); %#ok<AGROW>
    end
end
end

function runSel = load_run_selection(tsvFile, SUB, SES)
% LOAD_RUN_SELECTION  Load run_selection.tsv and return selections for this sub/ses.
%
% Returns struct with fields: T2w, task_<label>, reverse
% Each field is empty (no selection) or 'run-N'.
runSel = struct();

if isempty(tsvFile) || ~isfile(tsvFile)
    return;
end

fid = fopen(tsvFile, 'r');
if fid == -1, return; end

% Read and skip header
hdr = fgetl(fid);
if ~ischar(hdr), fclose(fid); return; end

while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    line = strtrim(line);
    if isempty(line), continue; end

    parts = strsplit(line, '\t');
    if numel(parts) < 5, continue; end

    fSub = strtrim(parts{1});
    fSes = strtrim(parts{2});
    fType = strtrim(parts{3});
    fSelected = strtrim(parts{5});

    if ~strcmp(fSub, SUB) || ~strcmp(fSes, SES)
        continue;
    end

    if isempty(fSelected), continue; end

    % Store selection: type -> selected run
    % Normalise field name for MATLAB struct
    fieldName = regexprep(fType, '[^a-zA-Z0-9]', '_');
    runSel.(fieldName) = fSelected;
end
fclose(fid);
end

function t2w = find_t2w(anatDir, SUB, SES, runSel)
% FIND_T2W  Find the T2w image, handling multi-run and .nii.gz.
%
% Priority:
%   1. If run_selection.tsv specifies a run, use that
%   2. If only one T2w file, use it
%   3. If multiple, use the last run (highest number) and warn

% Find all T2w NIfTI files
t2wFiles = all_matches(anatDir, ['^' regexptranslate('escape',SUB) '_' regexptranslate('escape',SES) '.*T2w\.nii']);

assert(~isempty(t2wFiles), 'No T2w file found in %s', anatDir);

if numel(t2wFiles) == 1
    t2w = t2wFiles{1};
    return;
end

% Multiple T2w files — check run selection
if isfield(runSel, 'T2w') && ~isempty(runSel.T2w)
    selRun = runSel.T2w;
    for i = 1:numel(t2wFiles)
        [~, fname] = fileparts(t2wFiles{i});
        % Handle .nii.gz double extension
        fname = regexprep(fname, '\.nii$', '');
        if contains(fname, selRun)
            t2w = t2wFiles{i};
            fprintf('  T2w: using selected %s (from run_selection.tsv)\n', selRun);
            return;
        end
    end
    warning('run_selection.tsv specifies %s for T2w but no matching file found. Using last run.', selRun);
end

% Default: use last run (sort and take last)
t2wFiles = sort(t2wFiles);
t2w = t2wFiles{end};
fprintf('  T2w: multiple runs found (%d), using last: %s\n', numel(t2wFiles), t2w);
end

function [roiL, roiR, hasRoiL, hasRoiR] = find_rois(anatDir, SUB, SES, CFG)
% FIND_ROIS  Find left and right ROI masks using flexible matching.
%
% Searches for files containing the ROI pattern as a substring,
% tolerating extra BIDS entities (run-N, acq-*, etc.).

roiL = ''; roiR = ''; hasRoiL = false; hasRoiR = false;

if ~CFG.DO_COREG_ROIS, return; end

patL = CFG.ROI_PATTERN_LEFT;
patR = CFG.ROI_PATTERN_RIGHT;

% Strategy 1: exact match  sub_ses + pattern
exact_L = fullfile(anatDir, [SUB '_' SES patL]);
exact_R = fullfile(anatDir, [SUB '_' SES patR]);

if isfile(exact_L)
    roiL = exact_L; hasRoiL = true;
else
    % Strategy 2: substring match (handles extra entities)
    escapedPat = regexptranslate('escape', patL);
    hit = first_match(anatDir, ['^' regexptranslate('escape',SUB) '_' regexptranslate('escape',SES) '.*' escapedPat]);
    if ~isempty(hit)
        roiL = hit; hasRoiL = true;
    end
end

if isfile(exact_R)
    roiR = exact_R; hasRoiR = true;
else
    escapedPat = regexptranslate('escape', patR);
    hit = first_match(anatDir, ['^' regexptranslate('escape',SUB) '_' regexptranslate('escape',SES) '.*' escapedPat]);
    if ~isempty(hit)
        roiR = hit; hasRoiR = true;
    end
end
end

function boldFile = find_bold(funcDir, SUB, SES, tLabel, runSel)
% FIND_BOLD  Find BOLD file for a task label, handling multi-run and .nii.gz.
%
% Looks for: sub-XX_ses-YY_task-<tLabel>[_run-N]_bold.nii[.gz]
% If multiple runs, uses run_selection or defaults to last run.

% Try exact match first (no run entity)
exact = fullfile(funcDir, sprintf('%s_%s_task-%s_bold.nii', SUB, SES, tLabel));
if isfile(exact)
    boldFile = exact;
    return;
end
% Try .nii.gz
exact_gz = [exact '.gz'];
if isfile(exact_gz)
    boldFile = exact_gz;
    return;
end

% Search with flexible regex: sub_ses_task-<label>[_anything]_bold.nii[.gz]
escLabel = regexptranslate('escape', tLabel);
pattern = ['^' regexptranslate('escape',SUB) '_' regexptranslate('escape',SES) '_task-' escLabel '(_[^/]*)?' '_bold\.nii'];
allBold = all_matches(funcDir, pattern);

if isempty(allBold)
    % Also try without strict anchoring (handles unexpected prefixes)
    pattern2 = ['task-' escLabel '.*_bold\.nii'];
    allBold = all_matches(funcDir, pattern2);
end

assert(~isempty(allBold), 'Cannot find BOLD for task-%s in %s', tLabel, funcDir);

if numel(allBold) == 1
    boldFile = allBold{1};
    return;
end

% Multiple runs found — check run selection
selKey = ['task_' tLabel];
if isfield(runSel, selKey) && ~isempty(runSel.(selKey))
    selRun = runSel.(selKey);
    for i = 1:numel(allBold)
        [~, fname] = fileparts(allBold{i});
        if contains(fname, selRun)
            boldFile = allBold{i};
            fprintf('  BOLD task-%s: using selected %s (from run_selection.tsv)\n', tLabel, selRun);
            return;
        end
    end
    warning('run_selection.tsv specifies %s for task-%s but no match found. Using last run.', selRun, tLabel);
end

% Default: sort and take last (highest run number)
allBold = sort(allBold);
boldFile = allBold{end};
fprintf('  BOLD task-%s: multiple runs (%d), using last: %s\n', tLabel, numel(allBold), boldFile);
end

function vols = expand_4d(niiFile)
V = spm_vol(niiFile);
vols = cell(numel(V),1);
for k=1:numel(V)
    vols{k} = sprintf('%s,%d', niiFile, k);
end
end

%

function vdm_out = calc_vdm_for_run_to_func(phasemap, magnitude1, epiRef, funcDir, fmapDir, CFG, taskN)
% Calculates a VDM matched to epiRef and ensures output ends up in funcDir
% as a deterministic name: vdm_task-<N>.nii
%
% Robust strategy:
% - Don't trust vdmflags.prefix (SPM/FieldMap sometimes ignores or overrides it)
% - Instead: snapshot existing vdm* files (func+fmap), run job, then diff to find newly created VDM.

% ---- Snapshot existing VDM candidates (before) ----
before = list_vdm_candidates(funcDir);
before = [before; list_vdm_candidates(fmapDir)];
before = unique(before);

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
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname      = sprintf('task-%d', taskN);
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
after = list_vdm_candidates(funcDir);
after = [after; list_vdm_candidates(fmapDir)];
after = unique(after);

newFiles = setdiff(after, before);

% If nothing appears in func/fmap, also check current dir (rare but happens)
if isempty(newFiles)
    newFiles = setdiff(list_vdm_candidates(pwd), before);
end

assert(~isempty(newFiles), 'VDM not found after FieldMap run for task-%d. Check where FieldMap writes outputs.', taskN);

% If multiple candidates were created, take the newest one
newest = pick_newest_file(newFiles);

% ---- Move/rename deterministically into funcDir ----
vdm_out = fullfile(funcDir, sprintf('vdm_task-%d.nii', taskN));

% If output is .img/.hdr pair, convert handling:
[~,~,ext] = fileparts(newest);
if strcmpi(ext,'.img')
    % Move both .img and .hdr, but keep .img name; SPM can read Analyze too.
    target_img = fullfile(funcDir, sprintf('vdm_task-%d.img', taskN));
    target_hdr = fullfile(funcDir, sprintf('vdm_task-%d.hdr', taskN));
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
function files = list_vdm_candidates(folder)
files = {};
if ~isfolder(folder), return; end
d = dir(folder);

for i = 1:numel(d)
    if d(i).isdir, continue; end
    n = d(i).name;

    % FieldMap typically uses vdm*.nii or vdm*.img/.hdr
    if ~isempty(regexp(n,'^vdm.*\.(nii|img)$','once'))
        files{end+1,1} = fullfile(folder,n); %#ok<AGROW>
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

% ===================== FSL TOPUP HELPER FUNCTIONS =====================

function reversePE_epi = find_reverse_pe_epi(subDir, SUB, SES, CFG, runSel)
% FIND_REVERSE_PE_EPI  Locate the reverse phase-encode EPI for topup.
%
% Searches fmap/ and/or func/ depending on TOPUP_REVERSE_PE_DIR setting.
% Uses flexible matching: the TOPUP_REVERSE_PE_PATTERN is treated as a
% substring that must appear somewhere in the NIfTI filename.  Extra BIDS
% entities (acq-*, run-*, etc.) between specifiers are tolerated.
% If multiple matches are found, run_selection.tsv is consulted.
%
% Search strategy per directory:
%   1) Exact BIDS name: <SUB>_<SES>_<pattern>.nii (with optional .gz)
%   2) Substring match: any .nii file whose name contains <pattern>
%   3) Broad fallback:  common reverse-PE names (dir-*_epi, task-reverse)

pattern = CFG.TOPUP_REVERSE_PE_PATTERN;
searchDir = CFG.TOPUP_REVERSE_PE_DIR;

fmapDir = fullfile(subDir, 'fmap');
funcDir = fullfile(subDir, 'func');

% Build ordered list of directories to search
switch searchDir
    case 'fmap'
        searchDirs = {fmapDir};
    case 'func'
        searchDirs = {funcDir};
    case 'auto'
        searchDirs = {};
        if isfolder(fmapDir), searchDirs{end+1} = fmapDir; end
        searchDirs{end+1} = funcDir;
    otherwise
        error('Invalid TOPUP_REVERSE_PE_DIR: %s (must be fmap, func, or auto)', searchDir);
end

% Collect all candidates across all search directories
allCandidates = {};

for di = 1:numel(searchDirs)
    d = searchDirs{di};
    if ~isfolder(d), continue; end

    % Strategy 1: exact BIDS name  <SUB>_<SES>_<pattern>.nii[.gz]
    candidate = fullfile(d, sprintf('%s_%s_%s.nii', SUB, SES, pattern));
    if isfile(candidate)
        allCandidates{end+1} = candidate; %#ok<AGROW>
    end
    candidate_gz = [candidate '.gz'];
    if isfile(candidate_gz)
        allCandidates{end+1} = candidate_gz; %#ok<AGROW>
    end

    % Strategy 2: substring match
    escapedPat = regexptranslate('escape', pattern);
    hits = all_matches(d, ['.*' escapedPat '.*\.nii']);
    for hi = 1:numel(hits)
        if ~any(strcmp(hits{hi}, allCandidates))
            allCandidates{end+1} = hits{hi}; %#ok<AGROW>
        end
    end
end

% Broad fallback if nothing found
if isempty(allCandidates)
    fallbackPatterns = { ...
        '.*dir-[A-Za-z]+_epi.*\.nii',    ...
        '.*task-reverse.*bold.*\.nii',    ...
        '.*task-reverse.*\.nii',          ...
        '.*_epi\.nii'                     ...
    };
    for di = 1:numel(searchDirs)
        d = searchDirs{di};
        if ~isfolder(d), continue; end
        for fi = 1:numel(fallbackPatterns)
            hits = all_matches(d, fallbackPatterns{fi});
            for hi = 1:numel(hits)
                if ~any(strcmp(hits{hi}, allCandidates))
                    allCandidates{end+1} = hits{hi}; %#ok<AGROW>
                    fprintf('  (matched via fallback pattern: %s)\n', fallbackPatterns{fi});
                end
            end
        end
    end
end

if isempty(allCandidates)
    searchedStr = strjoin(searchDirs, ', ');
    error(['Could not find reverse-PE EPI.\n' ...
           '  Pattern: %s\n' ...
           '  Searched: %s\n' ...
           '  TOPUP_REVERSE_PE_DIR: %s\n' ...
           '  Try adjusting TOPUP_REVERSE_PE_PATTERN in pipeline_config.cfg'], ...
           pattern, searchedStr, searchDir);
end

if numel(allCandidates) == 1
    reversePE_epi = allCandidates{1};
    return;
end

% Multiple candidates — check run_selection
if isfield(runSel, 'task_reverse') && ~isempty(runSel.task_reverse)
    selRun = runSel.task_reverse;
    for i = 1:numel(allCandidates)
        [~, fname] = fileparts(allCandidates{i});
        if contains(fname, selRun)
            reversePE_epi = allCandidates{i};
            fprintf('  Reverse-PE: using selected %s (from run_selection.tsv)\n', selRun);
            return;
        end
    end
    warning('run_selection.tsv specifies %s for reverse-PE but no match. Using first.', selRun);
end

% Default: use the first match
allCandidates = sort(allCandidates);
reversePE_epi = allCandidates{1};
fprintf('  Reverse-PE: multiple found (%d), using first: %s\n', numel(allCandidates), reversePE_epi);
end

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
fid = fopen(acqparamsFile, 'w');
assert(fid > 0, 'Cannot create acqparams file: %s', acqparamsFile);

% Convert PE direction strings to acqparams vectors
fwdVec = pe_dir_to_vector(CFG.TOPUP_PE_DIR_BOLD);
revVec = pe_dir_to_vector(CFG.TOPUP_PE_DIR_REVERSE);
readout = CFG.TOPUP_READOUT_SEC;

fprintf(fid, '%d %d %d %.6f\n', fwdVec(1), fwdVec(2), fwdVec(3), readout);
fprintf(fid, '%d %d %d %.6f\n', revVec(1), revVec(2), revVec(3), readout);
fclose(fid);
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

function outputFile = run_fsl_applytopup(boldFile, topupPrefix, outputFile, imainIndex, CFG)
% RUN_FSL_APPLYTOPUP  Apply topup distortion correction to a 4D BOLD file.
%
% Inputs:
%   boldFile    - Input 4D BOLD NIfTI
%   topupPrefix - Prefix from run_fsl_topup (contains fieldcoefs, movpar)
%   outputFile  - Desired output path (must end in .nii)
%   imainIndex  - Index into acqparams.txt for this image (1 = forward PE)
%   CFG         - Config struct with TOPUP_INTERP
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

acqparamsFile = fullfile(fileparts(topupPrefix), 'topup_acqparams.txt');

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

function vdm_out = convert_topup_to_vdm(topupPrefix, epiRef, funcDir, CFG, taskN)
% CONVERT_TOPUP_TO_VDM  Convert FSL topup field to SPM-compatible VDM.
%
% The topup field (in Hz) is converted to a voxel displacement map (in mm)
% using the formula: VDM = fieldmap_Hz * totalReadoutTime * voxelSize_PE
%
% The VDM is then written as a NIfTI in funcDir for use with SPM Realign&Unwarp.

fieldFile = sprintf('%s_field.nii.gz', topupPrefix);
assert(isfile(fieldFile), 'Topup field not found: %s', fieldFile);

% Convert to .nii for SPM (FSL defaults to .nii.gz)
fieldNii = fullfile(funcDir, sprintf('topup_field_task-%d.nii', taskN));
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
vdm_out = fullfile(funcDir, sprintf('vdm_task-%d.nii', taskN));
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

function niiFile = ensure_nii(inputFile)
% ENSURE_NII  Decompress .nii.gz to .nii if needed. Returns .nii path.
%
% SPM cannot read compressed NIfTI files. This function checks if the
% input is .nii.gz, decompresses it in-place using gunzip, and returns
% the path to the uncompressed .nii file.
%
% If the input is already .nii, it is returned unchanged.
% If a .nii already exists alongside the .nii.gz, the .nii is returned.

if isempty(inputFile)
    niiFile = inputFile;
    return;
end

% Check if it's .nii.gz
if endsWith(inputFile, '.nii.gz')
    niiFile = inputFile(1:end-3);  % strip .gz
    if isfile(niiFile)
        fprintf('  [ensure_nii] Already decompressed: %s\n', niiFile);
        return;
    end
    fprintf('  [ensure_nii] Decompressing: %s\n', inputFile);
    gunzip(inputFile);
    assert(isfile(niiFile), 'Decompression failed: %s not created', niiFile);
else
    niiFile = inputFile;
end
end

% ===================== BIDS JSON VALIDATION =====================

function validate_bids_params(funcDir, fmapDir, SUB, SES, CFG)
% VALIDATE_BIDS_PARAMS  Compare pipeline config against BIDS JSON sidecars.
%
% Reads JSON sidecar files from the BIDS dataset and compares key
% acquisition parameters against the values in pipeline_config.cfg.
% Prints warnings for any mismatches. Does NOT stop the pipeline —
% the config values are always used (in case the user intentionally
% overrides the BIDS values).

fprintf('\n--- Validating config against BIDS JSON sidecars ---\n');
nWarnings = 0;

% --- Find a task BOLD JSON to validate TR, PE direction, readout ---
tLabel = CFG.TASK_LABELS{1};
boldJsonPatterns = {
    fullfile(funcDir, sprintf('%s_%s_task-%s_bold.json', SUB, SES, tLabel))
    fullfile(funcDir, sprintf('%s_%s_task-%s_run-1_bold.json', SUB, SES, tLabel))
};
boldJson = '';
for i = 1:numel(boldJsonPatterns)
    if isfile(boldJsonPatterns{i})
        boldJson = boldJsonPatterns{i};
        break;
    end
end
% Fallback: any task BOLD JSON
if isempty(boldJson)
    hits = dir(fullfile(funcDir, sprintf('%s_%s_task-*_bold.json', SUB, SES)));
    hits = hits(~startsWith({hits.name}, '._'));
    if ~isempty(hits)
        boldJson = fullfile(funcDir, hits(1).name);
    end
end

% Initialize PE variables for use in both BOLD and reverse-PE validation
bidsPE = '';
bidsPE_hasSign = false;

if ~isempty(boldJson)
    bj = read_json(boldJson);
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
    bidsPE = '';
    bidsPE_hasSign = false;
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
            % Extract just the axis letter for axis-only comparison
            configAxis = regexprep(CFG.TOPUP_PE_DIR_BOLD, '-$', '');  % 'y-' -> 'y'
            if bidsPE_hasSign
                % Full direction available — exact match
                if ~strcmp(fslPE, CFG.TOPUP_PE_DIR_BOLD)
                    fprintf('  WARNING: PE direction mismatch — config TOPUP_PE_DIR_BOLD: %s, BIDS JSON: %s (=%s in FSL)\n', ...
                        CFG.TOPUP_PE_DIR_BOLD, bidsPE, fslPE);
                    nWarnings = nWarnings + 1;
                else
                    fprintf('  PE direction (BOLD): %s = %s in FSL (matches)\n', bidsPE, fslPE);
                end
            else
                % Only axis available — check axis matches, warn about sign
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

% --- Validate reverse-PE JSON for topup ---
if strcmp(CFG.PREPROC_MODE, 'topup')
    revPattern = CFG.TOPUP_REVERSE_PE_PATTERN;
    revJson = '';
    % Try common locations
    candidates = {};
    if isfolder(funcDir)
        hits = dir(fullfile(funcDir, ['*' revPattern '*.json']));
        for i = 1:numel(hits)
            candidates{end+1} = fullfile(funcDir, hits(i).name); %#ok<AGROW>
        end
    end
    if isfolder(fmapDir)
        hits = dir(fullfile(fmapDir, ['*' revPattern '*.json']));
        for i = 1:numel(hits)
            candidates{end+1} = fullfile(fmapDir, hits(i).name); %#ok<AGROW>
        end
    end
    % Also try standard naming
    stdCandidates = {
        fullfile(funcDir, sprintf('%s_%s_%s.json', SUB, SES, revPattern))
    };
    for i = 1:numel(stdCandidates)
        if isfile(stdCandidates{i})
            candidates{end+1} = stdCandidates{i}; %#ok<AGROW>
        end
    end
    if ~isempty(candidates)
        revJson = candidates{1};
    end

    if ~isempty(revJson)
        rj = read_json(revJson);
        fprintf('  Reverse-PE JSON: %s\n', revJson);

        % Check reverse PE direction
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
                % Full direction — exact match
                if ~isempty(fslRevPE) && ~strcmp(fslRevPE, CFG.TOPUP_PE_DIR_REVERSE)
                    fprintf('  WARNING: Reverse PE direction mismatch — config: %s, BIDS JSON: %s (=%s in FSL)\n', ...
                        CFG.TOPUP_PE_DIR_REVERSE, revPE, fslRevPE);
                    nWarnings = nWarnings + 1;
                else
                    fprintf('  PE direction (reverse): %s = %s in FSL (matches)\n', revPE, fslRevPE);
                end
            else
                % Axis only — check axis matches
                if ~strcmp(bidsRevAxis, configRevAxis)
                    fprintf('  WARNING: Reverse PE axis mismatch — config axis: %s, BIDS JSON: %s (=%s in FSL)\n', ...
                        configRevAxis, revPE, bidsRevAxis);
                    nWarnings = nWarnings + 1;
                else
                    fprintf('  PE axis (reverse): %s (matches). PhaseEncodingAxis only — cannot verify +/- sign.\n', bidsRevAxis);
                end
            end
        end

        % Validate PE vectors are opposite (only possible with signed directions)
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
        elseif ~isempty(bidsPE) && ~isempty(revPE) && ~(bidsPE_hasSign && revPE_hasSign)
            % Check at least that both are on the same axis (required for topup)
            fwdAxis = regexprep(bids_pe_to_fsl(bidsPE), '-$', '');
            revAxis = regexprep(bids_pe_to_fsl(revPE), '-$', '');
            if strcmp(fwdAxis, revAxis)
                fprintf('  Forward and reverse PE share axis %s (OK). Only PhaseEncodingAxis available — verify opposite signs manually.\n', fwdAxis);
            else
                fprintf('  WARNING: Forward PE axis (%s) and reverse PE axis (%s) are on DIFFERENT axes!\n', fwdAxis, revAxis);
                nWarnings = nWarnings + 1;
            end
        end
    else
        fprintf('  No reverse-PE JSON sidecar found — skipping reverse-PE validation.\n');
    end
end

% --- Validate fieldmap JSONs for realign_unwarp ---
if strcmp(CFG.PREPROC_MODE, 'realign_unwarp') && isfolder(fmapDir)
    phasediffJson = first_match(fmapDir, '.*phasediff.*\.json$');
    if ~isempty(phasediffJson)
        pj = read_json(phasediffJson);
        fprintf('  Phasediff JSON: %s\n', phasediffJson);

        % Check echo times
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
        fprintf('  No phasediff JSON found in fmap/ — skipping TE validation.\n');
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

function js = read_json(jsonFile)
% READ_JSON  Read a JSON file and return as MATLAB struct.
fid = fopen(jsonFile, 'r');
if fid == -1
    js = struct();
    return;
end
raw = fread(fid, inf, '*char')';
fclose(fid);
try
    js = jsondecode(raw);
catch
    js = struct();
    warning('Could not parse JSON: %s', jsonFile);
end
end

% ===================== PROVENANCE LOG =====================

function save_provenance_log(funcDir, SUB, SES, CFG, cfgFile, taskBoldFiles)
% SAVE_PROVENANCE_LOG  Write a JSON log of all preprocessing settings used.
%
% Saved as preproc_provenance.json in the func/ directory.
% Contains: preprocessing mode, all key parameters, timestamps,
% software versions, and input files.

logFile = fullfile(funcDir, 'preproc_provenance.json');
fprintf('\nSaving provenance log: %s\n', logFile);

prov = struct();

% --- Metadata ---
prov.created = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
prov.subject = SUB;
prov.session = SES;
prov.config_file = cfgFile;
prov.preproc_mode = CFG.PREPROC_MODE;

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

% --- Input files ---
prov.inputs.task_labels = CFG.TASK_LABELS;
prov.inputs.bold_files = taskBoldFiles;

% --- Acquisition parameters ---
if isfield(CFG, 'TR'), prov.acquisition.TR = CFG.TR; end
prov.acquisition.total_readout_ms = CFG.TOTAL_READOUT;
prov.acquisition.te_short_ms = CFG.TE_SHORT;
prov.acquisition.te_long_ms = CFG.TE_LONG;
prov.acquisition.blip_direction = CFG.BLIPDIR;
prov.acquisition.epi_based_fieldmap = CFG.EPI_BASED_FIELDMAP;

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
% Read the first BOLD JSON to store what BIDS says
tLabel = CFG.TASK_LABELS{1};
boldJsonCandidates = {
    fullfile(funcDir, sprintf('%s_%s_task-%s_bold.json', SUB, SES, tLabel))
    fullfile(funcDir, sprintf('%s_%s_task-%s_run-1_bold.json', SUB, SES, tLabel))
};
boldJson = '';
for i = 1:numel(boldJsonCandidates)
    if isfile(boldJsonCandidates{i})
        boldJson = boldJsonCandidates{i};
        break;
    end
end
if ~isempty(boldJson)
    bj = read_json(boldJson);
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
end
