function run_gridcat_analysis(cfgFile, subjectListFile)
% RUN_GRIDCAT_ANALYSIS  Run the GridCAT analysis (GLM1 + GLM2) for all subjects and sessions.
%
% Usage:
%   run_gridcat_analysis()                              % Auto-detects pipeline_config.cfg
%   run_gridcat_analysis('pipeline_config.cfg')        % Current directory
%   run_gridcat_analysis('/full/path/pipeline_config.cfg')  % Explicit path
%   run_gridcat_analysis(cfgFile, 'subses_list.txt')   % Only process subjects in list
%
% This function:
%   1. Reads all settings from pipeline_config.cfg (or auto-detects it).
%   2. Finds all EventData files, functional scans, motion regressors, and ROI masks.
%   3. Optionally filters to only process subject-sessions in subjectListFile.
%   4. Runs GLM1 (first-level GLM) for each subject-session.
%   5. Runs GLM2 (grid orientation analysis) using the ROI masks you specify.
%   6. Exports grid metrics (the main output you care about).
%
% Before running, configure the settings in pipeline_config.cfg.
% Everything else should work without changes.
%
% HOW TO RUN (on the HPC via SLURM):
%   Use the accompanying run_gridcat.sh script, or manually:
%     sbatch --time=08:00:00 --mem=32G --cpus-per-task=8 \
%            --wrap="module load matlab/R2024b; matlab -nodisplay -r run_gridcat_analysis"


%% =====================================================================
%%                      READ CONFIGURATION
%% =====================================================================

% If no config file is provided, auto-detect it in the same directory as this script
if nargin < 1 || isempty(cfgFile)
    scriptPath = mfilename('fullpath');
    scriptDir = fileparts(scriptPath);
    cfgFile = fullfile(scriptDir, 'pipeline_config.cfg');
end

% SubjectList filter: if provided, only process these subject-session pairs
if nargin < 2
    subjectListFile = '';
end

% Read the configuration from file
cfg = read_pipeline_config(cfgFile);

% Map config keys to local variables with defaults
% ---- Directories ----
NFS_OUTPUT_ROOT = cfg.OUTPUT_ROOT;
ROOT_DIR = fullfile(NFS_OUTPUT_ROOT, 'GLM_runauto');

% ---- Variant-aware output folder ----
% RUN_VARIANT in the config lets several GridCAT analyses coexist under
% the same OUTPUT_ROOT: empty -> GLM_output (default), else GLM_output_<variant>.
RUN_VARIANT_RAW = '';
if isfield(cfg, 'RUN_VARIANT') && ~isempty(cfg.RUN_VARIANT)
    if ischar(cfg.RUN_VARIANT) || isstring(cfg.RUN_VARIANT)
        RUN_VARIANT_RAW = char(cfg.RUN_VARIANT);
    else
        RUN_VARIANT_RAW = num2str(cfg.RUN_VARIANT);
    end
end
SAFE_VARIANT = regexprep(RUN_VARIANT_RAW, '[^A-Za-z0-9._-]+', '_');
SAFE_VARIANT = regexprep(SAFE_VARIANT, '^_+|_+$', '');
if isempty(SAFE_VARIANT)
    OUTPUT_FOLDER = 'GLM_output';
else
    OUTPUT_FOLDER = ['GLM_output_' SAFE_VARIANT];
end
OUTPUT_ROOT = fullfile(NFS_OUTPUT_ROOT, OUTPUT_FOLDER);

% ---- Local scratch staging (bypass NFS I/O bottleneck) ----
USE_LOCAL_SCRATCH = isfield(cfg, 'USE_LOCAL_SCRATCH') && strcmp(cfg.USE_LOCAL_SCRATCH, 'true');
LOCAL_SCRATCH_BASE = '';
if isfield(cfg, 'LOCAL_SCRATCH_BASE') && ~isempty(cfg.LOCAL_SCRATCH_BASE)
    LOCAL_SCRATCH_BASE = cfg.LOCAL_SCRATCH_BASE;
end

SPM_DIR = cfg.SPM_DIR;
GRIDCAT_DIR = cfg.GRIDCAT_DIR;

% CircStat directory: use from config, or default to fullfile(GRIDCAT_DIR, 'CircStat2012a')
if isfield(cfg, 'CIRCSTAT_DIR') && ~isempty(cfg.CIRCSTAT_DIR)
    CIRCSTAT_DIR = cfg.CIRCSTAT_DIR;
else
    CIRCSTAT_DIR = fullfile(GRIDCAT_DIR, 'CircStat2012a');
end

% ---- File prefixes ----
FUNC_PREFIX = cfg.FUNC_PREFIX;
if isempty(FUNC_PREFIX)
    FUNC_PREFIX = 'u';
end

ERC_PREFIX = cfg.ROI_PREFIX;
if isempty(ERC_PREFIX)
    ERC_PREFIX = 'r';
end

% ---- ROI mode ----
ROI_MODE = cfg.ROI_MODE;
if isempty(ROI_MODE)
    ROI_MODE = 'both';
end

% ---- Task filtering ----
INCLUDE_TASKS = {};  % Currently unused from config, can be extended

% ---- fMRI / GLM model settings ----
GLM_TR = cfg.TR;
X_FOLD = cfg.X_FOLD_SYMMETRY;
MASKING_THR = cfg.MASKING_THRESHOLD;
MICRO_ONSET = cfg.MICROTIME_ONSET;
MICRO_RES = cfg.MICROTIME_RESOLUTION;
DERIVATIVES = cfg.DERIVATIVES;  % Expected as [0 0] from config
DISP_DESIGN = cfg.DISPLAY_DESIGN;
HPF_DEFAULT = cfg.HPF_CUTOFF;

% ---- GridCAT-specific settings ----
eventUsage_GLM1 = cfg.EVENT_USAGE_GLM1;
eventUsage_GLM2 = cfg.EVENT_USAGE_GLM2;
keepUnusedGridEvents = cfg.KEEP_UNUSED_GRID_EVENTS;
useVoxelWeighting = cfg.USE_VOXEL_WEIGHTING;
avgOriAcrossRuns_flag = cfg.AVG_ORI_ACROSS_RUNS;
GLM2_gridRegressorMethod = cfg.GLM2_REGRESSOR_METHOD;

% ---- Parallel processing ----
MAX_WORKERS = cfg.MAX_WORKERS;

% ---- Error handling (convert string to boolean) ----
FAIL_ON_MISSING = strcmp(cfg.FAIL_ON_MISSING, 'true');


%% =====================================================================
%%              PRINT CONFIGURATION SUMMARY
%% =====================================================================

fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('  GRIDCAT ANALYSIS - CONFIGURATION SUMMARY\n');
fprintf('%s\n', repmat('=', 1, 70));
fprintf('Config file: %s\n\n', cfgFile);

fprintf('DIRECTORIES:\n');
if isempty(RUN_VARIANT_RAW)
    fprintf('  Run variant:        <none>\n');
else
    fprintf('  Run variant:        %s\n', RUN_VARIANT_RAW);
end
fprintf('  Root directory:     %s\n', ROOT_DIR);
fprintf('  Output root:        %s\n', OUTPUT_ROOT);
fprintf('  SPM directory:      %s\n', SPM_DIR);
fprintf('  GridCAT directory:  %s\n', GRIDCAT_DIR);
fprintf('  CircStat directory: %s\n\n', CIRCSTAT_DIR);

fprintf('FILE PREFIXES:\n');
fprintf('  Functional prefix:  %s\n', FUNC_PREFIX);
fprintf('  ROI prefix:         %s\n\n', ERC_PREFIX);

fprintf('ROI & ANALYSIS SETTINGS:\n');
fprintf('  ROI mode:           %s\n', ROI_MODE);
fprintf('  TR:                 %.2f seconds\n', GLM_TR);
fprintf('  X-fold symmetry:    %d\n', X_FOLD);
fprintf('  Masking threshold:  %.2f\n', MASKING_THR);
fprintf('  Microtime onset:    %d\n', MICRO_ONSET);
fprintf('  Microtime resolution: %d\n', MICRO_RES);
fprintf('  Derivatives:        [%s]\n', sprintf('%d ', DERIVATIVES));
fprintf('  Display design:     %d\n', DISP_DESIGN);
fprintf('  HPF cutoff:         %d seconds\n\n', HPF_DEFAULT);

fprintf('GRIDCAT SETTINGS:\n');
fprintf('  Event usage GLM1:   %d\n', eventUsage_GLM1);
fprintf('  Event usage GLM2:   %d\n', eventUsage_GLM2);
fprintf('  Keep unused events: %d\n', keepUnusedGridEvents);
fprintf('  Use voxel weighting: %d\n', useVoxelWeighting);
fprintf('  Avg orientation across runs: %d\n', avgOriAcrossRuns_flag);
fprintf('  GLM2 regressor method: %s\n\n', GLM2_gridRegressorMethod);

fprintf('PARALLEL PROCESSING:\n');
fprintf('  Max workers:        %d\n', MAX_WORKERS);
fprintf('  Fail on missing:    %s\n', mat2str(FAIL_ON_MISSING));
fprintf('%s\n\n', repmat('=', 1, 70));

% Prompt user to continue if running interactively
if usejava('desktop')
    fprintf('Press Enter to continue or Ctrl+C to abort...\n');
    input('');
end


%% =====================================================================
%%                      SETUP (no need to edit below)
%% =====================================================================

% Check that all paths exist
assert(exist(ROOT_DIR, 'dir') == 7,     'ROOT_DIR not found: %s', ROOT_DIR);
assert(exist(SPM_DIR, 'dir') == 7,      'SPM_DIR not found: %s', SPM_DIR);
assert(exist(GRIDCAT_DIR, 'dir') == 7,  'GRIDCAT_DIR not found: %s', GRIDCAT_DIR);
assert(exist(CIRCSTAT_DIR, 'dir') == 7, 'CIRCSTAT_DIR not found: %s', CIRCSTAT_DIR);

if ~exist(OUTPUT_ROOT, 'dir'); mkdir(OUTPUT_ROOT); end

% Add toolboxes to MATLAB path
addpath(SPM_DIR);
addpath(GRIDCAT_DIR);
addpath(CIRCSTAT_DIR);
rehash toolboxcache;

% Start SPM in command-line mode (no GUI)
spm('defaults', 'fmri');
spm_jobman('initcfg');
spm_get_defaults('cmdline', true);

% Print which toolbox versions are being used
fprintf('SPM:      %s\n', which('spm'));
fprintf('GridCAT:  %s\n', which('specifyGLM'));
fprintf('CircStat: %s\n', which('circ_rtest'));

% Validate ROI_MODE
validModes = {'bilat_only', 'lr_only', 'both'};
assert(any(strcmpi(ROI_MODE, validModes)), ...
    'ROI_MODE must be one of: %s (you set: %s)', strjoin(validModes, ', '), ROI_MODE);

fprintf('\nROI_MODE    = %s\n', ROI_MODE);
fprintf('ROOT_DIR    = %s\n', ROOT_DIR);
fprintf('OUTPUT_ROOT = %s\n\n', OUTPUT_ROOT);


%% =====================================================================
%%              STAGE DATA TO LOCAL SCRATCH (optional)
%% =====================================================================
% When USE_LOCAL_SCRATCH=true, copies input data from NFS to node-local
% storage before processing. This eliminates the NFS I/O bottleneck for
% the many small .nii read/writes that SPM GLM estimation produces.
% Results are copied back to NFS at the end.

NFS_OUTPUT_DIR = OUTPUT_ROOT;  % Remember NFS destination for results
localScratchDir = '';

if USE_LOCAL_SCRATCH
    % Determine local scratch base: config, SLURM TMPDIR, or /tmp
    if ~isempty(LOCAL_SCRATCH_BASE)
        scratchBase = LOCAL_SCRATCH_BASE;
    elseif ~isempty(getenv('TMPDIR'))
        scratchBase = getenv('TMPDIR');
    else
        scratchBase = '/tmp';
    end

    % Create unique scratch directory for this job
    jobId = getenv('SLURM_JOB_ID');
    if isempty(jobId); jobId = sprintf('%d', feature('getpid')); end
    localScratchDir = fullfile(scratchBase, sprintf('gridcat_%s', jobId));

    localInput  = fullfile(localScratchDir, 'GLM_runauto');
    localOutput = fullfile(localScratchDir, OUTPUT_FOLDER);

    fprintf('=== LOCAL SCRATCH STAGING ===\n');
    fprintf('  Source (NFS): %s\n', ROOT_DIR);
    fprintf('  Local copy:   %s\n', localInput);
    fprintf('  Local output: %s\n', localOutput);

    % Diagnostic: show what filesystem /tmp actually is
    [~, fsInfo] = system('df -Th /tmp 2>/dev/null | tail -1');
    fprintf('  /tmp filesystem: %s', strtrim(fsInfo));
    fprintf('\n  TMPDIR env:      "%s"\n', getenv('TMPDIR'));
    fprintf('  SLURM_TMPDIR:    "%s"\n', getenv('SLURM_TMPDIR'));
    fprintf('  Scratch base:    %s\n', scratchBase);

    % Check if /tmp is actually local (not NFS) — skip staging if it's NFS
    [~, fsType] = system('stat -f -c %T /tmp 2>/dev/null');
    fsType = strtrim(fsType);
    if contains(fsType, 'nfs', 'IgnoreCase', true)
        warning('/tmp is NFS-backed (%s). Local scratch staging would not help. Falling back to NFS.', fsType);
        USE_LOCAL_SCRATCH = false;
    end

    if USE_LOCAL_SCRATCH
        tStage = tic;

        % Create the scratch directory tree BEFORE rsync
        [rc, msg] = system(sprintf('mkdir -p "%s" "%s"', localInput, localOutput));
        if rc ~= 0
            warning('Cannot create scratch directories (rc=%d): %s\nFalling back to NFS.', rc, msg);
            USE_LOCAL_SCRATCH = false;
        end
    end

    if USE_LOCAL_SCRATCH
        % Copy input data to local storage using tar pipe (faster than rsync
        % for many small files over NFS — single streaming read instead of
        % 30k+ individual file open/stat/read/close operations)
        tarCmd = sprintf('tar -C "%s" -cf - . | tar -C "%s" -xf -', ROOT_DIR, localInput);
        fprintf('  Copying input data to local scratch (tar pipe)...\n');
        [rc, msg] = system(tarCmd);
        if rc ~= 0
            warning('Copy to scratch failed (rc=%d): %s\nFalling back to NFS.', rc, msg);
            USE_LOCAL_SCRATCH = false;
        else
            fprintf('  Staged %.1f GB in %.1f seconds.\n', ...
                dir_size_gb(localInput), toc(tStage));

            % Redirect paths to local scratch
            ROOT_DIR    = localInput;
            OUTPUT_ROOT = localOutput;

            fprintf('  ROOT_DIR    -> %s\n', ROOT_DIR);
            fprintf('  OUTPUT_ROOT -> %s\n', OUTPUT_ROOT);
        end
    end
    fprintf('\n');
end


%% =====================================================================
%%                    FIND ALL FILES
%% =====================================================================
% Scan each known subdirectory ONCE upfront and index files in memory.
% This avoids repeated recursive dir('**',...) calls over 30k+ files.

fprintf('Indexing files...\n');
tIndex = tic;

% --- Load subject list filter early (before file discovery) ---
allowedPairs = containers.Map('KeyType', 'char', 'ValueType', 'logical');
hasSubjectFilter = false;
if ~isempty(subjectListFile)
    if ~isfile(subjectListFile)
        error('SubjectList file not found: %s', subjectListFile);
    end
    allowedPairs = load_allowed_pairs(subjectListFile);
    hasSubjectFilter = true;
end

% --- Scan each subdirectory once (non-recursive, fast) ---
eventDir = fullfile(ROOT_DIR, 'EventFiles', 'Event_tables');
rpDir    = fullfile(ROOT_DIR, 'rp_txt');
funcDir  = fullfile(ROOT_DIR, 'functional_scans_split');
roiDir   = fullfile(ROOT_DIR, 'ROI');

% Fall back to recursive scan only if the known subdirectories don't exist
if exist(eventDir, 'dir') && exist(rpDir, 'dir') && exist(funcDir, 'dir')
    eventFiles = dir(fullfile(eventDir, '*EventData.txt'));
    rpAll      = dir(fullfile(rpDir, 'rp*.txt'));
    funcAll    = dir(fullfile(funcDir, [FUNC_PREFIX '*_bold_*.nii']));
    roiAll     = dir(fullfile(roiDir, [ERC_PREFIX '*.nii']));
    fprintf('  Fast index: %d events, %d rp, %d func, %d ROI files (%.1fs)\n', ...
        numel(eventFiles), numel(rpAll), numel(funcAll), numel(roiAll), toc(tIndex));
else
    warning('Known subdirectory layout not found — falling back to recursive scan.');
    eventFiles = dir(fullfile(ROOT_DIR, '**', '*EventData.txt'));
    rpAll      = dir(fullfile(ROOT_DIR, '**', 'rp*.txt'));
    funcAll    = dir(fullfile(ROOT_DIR, '**', [FUNC_PREFIX '*_bold_*.nii']));
    roiAll     = dir(fullfile(ROOT_DIR, '**', [ERC_PREFIX '*.nii']));
    fprintf('  Recursive index: %d events, %d rp, %d func, %d ROI files (%.1fs)\n', ...
        numel(eventFiles), numel(rpAll), numel(funcAll), numel(roiAll), toc(tIndex));
end

% Remove macOS resource fork files
eventFiles = eventFiles(~startsWith({eventFiles.name}, '._'));
rpAll      = rpAll(~startsWith({rpAll.name}, '._'));
funcAll    = funcAll(~startsWith({funcAll.name}, '._'));
roiAll     = roiAll(~startsWith({roiAll.name}, '._'));

if isempty(eventFiles)
    error(['No *EventData.txt files found under ROOT_DIR.\n' ...
           'Expected location: %s'], eventDir);
end
fprintf('Found %d EventData files.\n', numel(eventFiles));

% --- Build lookup maps for rp and functional files (indexed by name) ---
rpNames   = {rpAll.name};
rpFolders = {rpAll.folder};
funcNames   = {funcAll.name};
funcFolders = {funcAll.folder};

% --- Which run of each task was preprocessed --------------------------------
% prepare_gridcat_directory.m flattens one run per task into GLM_runauto and
% records which one in run_manifest.tsv. Event tables are matched against that:
% a table naming a run the pipeline did not preprocess is not the right table
% for this data, and a table naming no run at all belongs to whichever run was.
runManifest = load_run_manifest(fullfile(ROOT_DIR, 'run_manifest.tsv'));
if runManifest.Count > 0
    fprintf('Run manifest: %d task entries.\n', runManifest.Count);
else
    fprintf(['No run_manifest.tsv in %s — event tables will be matched on\n' ...
             '  subject/session/task alone. Re-run the GridCAT preparation to\n' ...
             '  get one.\n'], ROOT_DIR);
end

% --- Parse each EventData file to build a list of runs ---
% The run entity is optional and captured separately: event tables are often
% exported without one.
patEvent = ['^(?<sub>sub-[^_]+)_(?<ses>ses-[^_]+)_(?<task>task-[A-Za-z0-9]+)' ...
            '(?<rest>.*)EventData\.txt$'];

runs = struct('sub', {}, 'ses', {}, 'task', {}, ...
              'eventFile', {}, 'addRegFile', {}, 'functionalScans', {});
claimed = containers.Map('KeyType', 'char', 'ValueType', 'char');

for i = 1:numel(eventFiles)
    tok = regexp(eventFiles(i).name, patEvent, 'names');
    if isempty(tok); continue; end

    % Apply task filter
    if ~isempty(INCLUDE_TASKS) && ~any(strcmp(tok.task, INCLUDE_TASKS))
        continue;
    end

    % Apply subject list filter early (skip files we won't need)
    if hasSubjectFilter
        pairKey = [tok.sub ' ' tok.ses];
        if ~isKey(allowedPairs, pairKey)
            continue;
        end
    end

    taskKey  = sprintf('%s|%s|%s', tok.sub, tok.ses, tok.task);
    eventRun = run_entity_of(tok.rest);

    % -- Does this event table describe the run that was preprocessed? --
    if isKey(runManifest, taskKey)
        usedRun = runManifest(taskKey);
        if ~isempty(eventRun) && ~strcmp(eventRun, usedRun)
            fprintf(['  Skipping %s: it is the event table for %s, but %s was\n' ...
                     '    preprocessed for %s %s %s.\n'], ...
                    eventFiles(i).name, eventRun, usedRun, tok.sub, tok.ses, tok.task);
            continue;
        end
    end

    % -- Two tables for one task cannot both be right --
    if isKey(claimed, taskKey)
        msg = sprintf([ ...
            'Two event tables claim %s %s %s:\n  %s\n  %s\n' ...
            'Only one run per task is preprocessed, so only one of these applies.\n' ...
            'Put the run entity in the filename (..._run-2_EventData.txt) so they\n' ...
            'can be told apart, or remove the one that does not belong.'], ...
            tok.sub, tok.ses, tok.task, claimed(taskKey), eventFiles(i).name);
        if FAIL_ON_MISSING; error(msg); else; warning(msg); end
        continue;
    end

    run.sub  = tok.sub;
    run.ses  = tok.ses;
    run.task = tok.task;
    run.eventFile = fullfile(eventFiles(i).folder, eventFiles(i).name);

    % --- Find motion regressor from pre-indexed list ---
    % The preparation writes rp_<sub>_<ses>_<task>.txt, so try that exact name
    % before falling back to a substring search over older layouts.
    rpWanted = sprintf('rp_%s_%s_%s.txt', tok.sub, tok.ses, tok.task);
    rpIdx = find(strcmp(rpNames, rpWanted), 1);
    if isempty(rpIdx)
        rpToken = sprintf('%s_%s_%s_', tok.sub, tok.ses, tok.task);
        rpMatch = contains(rpNames, rpToken) | ...
                  contains(rpNames, sprintf('%s_%s_%s.', tok.sub, tok.ses, tok.task));
        rpIdx = find(rpMatch);
    end

    if isempty(rpIdx)
        msg = sprintf('No motion regressor found for %s %s %s', tok.sub, tok.ses, tok.task);
        if FAIL_ON_MISSING; error(msg); else; warning(msg); continue; end
    end
    rpFull = fullfile(rpFolders(rpIdx), rpNames(rpIdx));
    [~, sortIdx] = sort(cellfun(@numel, rpFull));
    run.addRegFile = rpFull{sortIdx(1)};

    % --- Find 3D functional volumes from pre-indexed list ---
    funcToken = sprintf('%s%s_%s_%s_bold', FUNC_PREFIX, tok.sub, tok.ses, tok.task);
    funcMatch = startsWith(funcNames, [funcToken '_']);

    if ~any(funcMatch)
        msg = sprintf('No functional scans found for %s %s %s (prefix: %s)', ...
                       tok.sub, tok.ses, tok.task, funcToken);
        if FAIL_ON_MISSING; error(msg); else; warning(msg); continue; end
    end
    funcIdx = find(funcMatch);
    funcFull = fullfile(funcFolders(funcIdx), funcNames(funcIdx));
    funcFull = sort_by_trailing_number(funcFull);
    run.functionalScans = funcFull(:);

    claimed(taskKey) = eventFiles(i).name;
    runs(end+1) = run; %#ok<AGROW>
end

if isempty(runs)
    error('EventData files were found, but none matched the naming pattern or task filter.');
end

fprintf('Matched %d runs after filtering.\n', numel(runs));

% --- Group runs by subject-session ---
keys = strcat({runs.sub}, '|', {runs.ses});
[uKeys, ~, keyIdx] = unique(keys, 'stable');
fprintf('Found %d runs across %d subject-session pairs.\n\n', numel(runs), numel(uKeys));

% --- Find ROI masks from pre-indexed list ---
roiNames   = {roiAll.name};
roiFolders = {roiAll.folder};

sessInfo = repmat(struct('sub', '', 'ses', '', 'roi', ...
    struct('bilat', '', 'left', '', 'right', '')), numel(uKeys), 1);

for k = 1:numel(uKeys)
    parts = strsplit(uKeys{k}, '|');
    sub = parts{1}; ses = parts{2};
    sessInfo(k).sub = sub;
    sessInfo(k).ses = ses;

    if strcmpi(ROI_MODE, 'bilat_only') || strcmpi(ROI_MODE, 'both')
        sessInfo(k).roi.bilat = find_roi_from_index(roiNames, roiFolders, sub, ses, 'bilat', FAIL_ON_MISSING);
    end
    if strcmpi(ROI_MODE, 'lr_only') || strcmpi(ROI_MODE, 'both')
        sessInfo(k).roi.left  = find_roi_from_index(roiNames, roiFolders, sub, ses, 'left', FAIL_ON_MISSING);
        sessInfo(k).roi.right = find_roi_from_index(roiNames, roiFolders, sub, ses, 'right', FAIL_ON_MISSING);
    end
end

fprintf('File indexing complete (%.1f seconds).\n\n', toc(tIndex));


%% =====================================================================
%%                    SET UP PARALLEL PROCESSING
%% =====================================================================

useParallel = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));

if useParallel
    % Cap workers at number of sessions (no benefit from idle workers)
    nWorkers = min(MAX_WORKERS, numel(uKeys));

    % Per-job storage location so concurrent GridCAT runs do not collide.
    % SLURM_JOB_ID is unique per submission; fall back to PID for interactive runs.
    poolJobId = getenv('SLURM_JOB_ID');
    if isempty(poolJobId); poolJobId = sprintf('pid%d', feature('getpid')); end
    jobStorage = fullfile(getenv('HOME'), 'matlab_job_storage', poolJobId);
    if ~exist(jobStorage, 'dir'); mkdir(jobStorage); end
    % Clean the per-job storage on exit (best effort).
    cleanupJobStorage = onCleanup(@() cleanup_dir_safe(jobStorage));

    pc = parcluster('local');
    pc.JobStorageLocation = jobStorage;
    pc.NumWorkers = nWorkers;

    pool = gcp('nocreate');
    if ~isempty(pool); delete(pool); end

    fprintf('Starting parallel pool with %d workers (for %d sessions)...\n', nWorkers, numel(uKeys));
    pool = parpool(pc, nWorkers);

    % Initialize toolboxes on all workers using parfevalOnAll (non-blocking)
    spmDir_  = char(SPM_DIR);
    gridDir_ = char(GRIDCAT_DIR);
    circDir_ = char(CIRCSTAT_DIR);

    cmd = sprintf([ ...
        'addpath(''%s''); addpath(''%s''); addpath(''%s''); ' ...
        'spm(''defaults'',''fmri''); spm_jobman(''initcfg''); ' ...
        'spm_get_defaults(''cmdline'', true); '], ...
        strrep(spmDir_, '''', ''''''), ...
        strrep(gridDir_, '''', ''''''), ...
        strrep(circDir_, '''', ''''''));
    pctRunOnAll(cmd);

    fprintf('Parallel pool ready (%d workers).\n\n', nWorkers);
else
    fprintf('Parallel Computing Toolbox not available. Running one session at a time.\n\n');
end


%% =====================================================================
%%                    RUN THE ANALYSIS
%% =====================================================================

status = repmat(struct('sub', '', 'ses', '', 'ok', false, 'err', ''), numel(uKeys), 1);
tStart = tic;

if useParallel
    parfor k = 1:numel(uKeys)
        try
            thisRuns = runs(keyIdx == k);
            thisROI  = sessInfo(k).roi;
            run_one_session(thisRuns, thisROI, OUTPUT_ROOT, ROI_MODE, ...
                GLM_TR, X_FOLD, MASKING_THR, MICRO_ONSET, MICRO_RES, ...
                DERIVATIVES, DISP_DESIGN, HPF_DEFAULT, ...
                eventUsage_GLM1, eventUsage_GLM2, keepUnusedGridEvents, ...
                useVoxelWeighting, avgOriAcrossRuns_flag, GLM2_gridRegressorMethod);
            status(k) = struct('sub', thisRuns(1).sub, 'ses', thisRuns(1).ses, 'ok', true, 'err', '');
        catch ME
            status(k) = struct('sub', sessInfo(k).sub, 'ses', sessInfo(k).ses, ...
                'ok', false, 'err', getReport(ME, 'extended', 'hyperlinks', 'off'));
            write_error_log(OUTPUT_ROOT, sessInfo(k).sub, sessInfo(k).ses, status(k).err);
        end
    end
else
    for k = 1:numel(uKeys)
        try
            thisRuns = runs(keyIdx == k);
            thisROI  = sessInfo(k).roi;
            run_one_session(thisRuns, thisROI, OUTPUT_ROOT, ROI_MODE, ...
                GLM_TR, X_FOLD, MASKING_THR, MICRO_ONSET, MICRO_RES, ...
                DERIVATIVES, DISP_DESIGN, HPF_DEFAULT, ...
                eventUsage_GLM1, eventUsage_GLM2, keepUnusedGridEvents, ...
                useVoxelWeighting, avgOriAcrossRuns_flag, GLM2_gridRegressorMethod);
            status(k) = struct('sub', thisRuns(1).sub, 'ses', thisRuns(1).ses, 'ok', true, 'err', '');
        catch ME
            status(k) = struct('sub', sessInfo(k).sub, 'ses', sessInfo(k).ses, ...
                'ok', false, 'err', getReport(ME, 'extended', 'hyperlinks', 'off'));
            write_error_log(OUTPUT_ROOT, sessInfo(k).sub, sessInfo(k).ses, status(k).err);
        end
    end
end

elapsed = toc(tStart);


%% =====================================================================
%%              COPY RESULTS BACK FROM LOCAL SCRATCH
%% =====================================================================

if USE_LOCAL_SCRATCH && ~isempty(localScratchDir)
    fprintf('\n=== COPYING RESULTS BACK TO NFS ===\n');
    fprintf('  Local output: %s\n', OUTPUT_ROOT);
    fprintf('  NFS target:   %s\n', NFS_OUTPUT_DIR);

    tCopy = tic;
    % Ensure NFS output directory exists
    if ~exist(NFS_OUTPUT_DIR, 'dir'); mkdir(NFS_OUTPUT_DIR); end
    copyBack = sprintf('tar -C "%s" -cf - . | tar -C "%s" -xf -', OUTPUT_ROOT, NFS_OUTPUT_DIR);
    [rc, msg] = system(copyBack);
    if rc ~= 0
        warning('Copy-back failed (rc=%d): %s\nResults remain in: %s', rc, msg, OUTPUT_ROOT);
    else
        fprintf('  Copied results back in %.1f seconds.\n', toc(tCopy));

        % Clean up local scratch
        fprintf('  Cleaning up local scratch: %s\n', localScratchDir);
        rmdir(localScratchDir, 's');
    end

    % Point output path back to NFS for the summary
    OUTPUT_ROOT = NFS_OUTPUT_DIR;
end


%% =====================================================================
%%                    PRINT RESULTS SUMMARY
%% =====================================================================

nOK  = sum([status.ok]);
nBad = numel(status) - nOK;

fprintf('\n===========================\n');
fprintf('  ANALYSIS COMPLETE\n');
fprintf('===========================\n');
fprintf('Time:   %.1f minutes\n', elapsed / 60);
fprintf('OK:     %d / %d sessions\n', nOK, numel(status));
fprintf('Failed: %d / %d sessions\n', nBad, numel(status));
fprintf('Output: %s\n', OUTPUT_ROOT);

if nBad > 0
    fprintf('\nFailed sessions (check _logs folder for details):\n');
    for i = 1:numel(status)
        if ~status(i).ok
            fprintf('  - %s %s\n', status(i).sub, status(i).ses);
        end
    end
end

write_status_csv(OUTPUT_ROOT, status);
fprintf('\nDone.\n');

end


%% =====================================================================
%%           HELPER FUNCTIONS (you should not need to edit these)
%% =====================================================================

function run_one_session(thisRuns, roi, OUTPUT_ROOT, ROI_MODE, ...
    GLM_TR, X_FOLD, MASKING_THR, MICRO_ONSET, MICRO_RES, ...
    DERIVATIVES, DISP_DESIGN, HPF_DEFAULT, ...
    eventUsage_GLM1, eventUsage_GLM2, keepUnusedGridEvents, ...
    useVoxelWeighting, avgOriAcrossRuns_flag, GLM2_gridRegressorMethod)
% Processes one subject-session: runs GLM1, then GLM2 for each ROI.

    sub = thisRuns(1).sub;
    ses = thisRuns(1).ses;
    nRuns = numel(thisRuns);

    fprintf('\n=== %s %s (%d runs) ===\n', sub, ses, nRuns);

    % Build the cfg structure that GridCAT expects
    cfg = struct();
    for r = 1:nRuns
        cfg.rawData.run(r).functionalScans           = thisRuns(r).functionalScans;
        cfg.rawData.run(r).eventTable_file            = thisRuns(r).eventFile;
        cfg.rawData.run(r).additionalRegressors_file  = thisRuns(r).addRegFile;
    end

    cfg.GLM.TR                  = GLM_TR;
    cfg.GLM.xFoldSymmetry       = X_FOLD;
    cfg.GLM.maskingThreshold    = MASKING_THR;
    cfg.GLM.microtimeOnset      = MICRO_ONSET;
    cfg.GLM.microtimeResolution = MICRO_RES;
    cfg.GLM.HPF_perRun          = HPF_DEFAULT * ones(1, nRuns);
    cfg.GLM.derivatives         = DERIVATIVES;
    cfg.GLM.dispDesignMatrix    = DISP_DESIGN;
    cfg.GLM.keepUnusedGridEvents = keepUnusedGridEvents;

    % --- Step 1: GLM1 ---
    GLM1dir = fullfile(OUTPUT_ROOT, sprintf('%s_%s_GLM1', sub, ses));
    cfg.GLMnr = 1;
    cfg.GLM.dataDir = GLM1dir;
    cfg.GLM.eventUsageSpecifier = eventUsage_GLM1;

    fprintf('  GLM1 -> %s\n', GLM1dir);
    specifyGLM(cfg);
    estimateGLM(cfg);

    % --- Step 2: GLM2 for each ROI ---
    roiList = build_roi_list(roi, ROI_MODE, sub, ses);

    for j = 1:numel(roiList)
        roiLabel = roiList{j}.label;
        roiPath  = roiList{j}.path;

        cfg2 = cfg;
        cfg2.GLMnr = 2;
        cfg2.GLM.eventUsageSpecifier = eventUsage_GLM2;
        cfg2.GLM.GLM1_resultsDir = GLM1dir;
        cfg2.GLM.GLM2_roiMask_calcMeanGridOri      = {roiPath};
        cfg2.GLM.GLM2_useWeightingForVoxels         = useVoxelWeighting;
        cfg2.GLM.GLM2_averageMeanGridOriAcrossRuns   = avgOriAcrossRuns_flag;
        cfg2.GLM.GLM2_gridRegressorMethod            = GLM2_gridRegressorMethod;

        GLM2dir = fullfile(OUTPUT_ROOT, sprintf('%s_%s_GLM2_%s', sub, ses, roiLabel));
        cfg2.GLM.dataDir = GLM2dir;

        fprintf('  GLM2 (%s) -> %s\n', roiLabel, GLM2dir);
        specifyGLM(cfg2);
        estimateGLM(cfg2);

        % Export the grid metrics (this is the main result)
        outFile = {fullfile(OUTPUT_ROOT, sprintf('%s_%s_%s_grid_metrics.txt', sub, ses, roiLabel))};
        gridMetric_export({roiPath}, GLM1dir, GLM2dir, outFile);
    end

    fprintf('  Done %s %s.\n', sub, ses);
end


function roiList = build_roi_list(roi, ROI_MODE, sub, ses)
% Builds a list of ROI masks to use for GLM2, based on ROI_MODE.
    roiList = {};

    if strcmpi(ROI_MODE, 'bilat_only') || strcmpi(ROI_MODE, 'both')
        if ~isempty(roi.bilat)
            roiList{end+1} = struct('label', 'ErC-bilat', 'path', roi.bilat);
        end
    end

    if strcmpi(ROI_MODE, 'lr_only') || strcmpi(ROI_MODE, 'both')
        if ~isempty(roi.left)
            roiList{end+1} = struct('label', 'ErC-left', 'path', roi.left);
        end
        if ~isempty(roi.right)
            roiList{end+1} = struct('label', 'ErC-right', 'path', roi.right);
        end
    end

    if isempty(roiList)
        error('No ROI masks found for %s %s with ROI_MODE=%s', sub, ses, ROI_MODE);
    end
end


function roiPath = find_roi_from_index(roiNames, roiFolders, sub, ses, side, failFast)
% Finds a specific ROI mask from the pre-indexed ROI file list.
% Filters for names containing the subject, session, 'ErC', and the side label.
    roiPath = '';

    if ~ismember(lower(side), {'bilat', 'left', 'right'})
        error('Unknown ROI side: %s (expected bilat, left, or right)', side);
    end

    % Filter from pre-indexed list: must contain sub, ses, ErC, and side
    hasSub  = contains(roiNames, sub);
    hasSes  = contains(roiNames, ses);
    isErC   = ~cellfun(@isempty, regexpi(roiNames, 'ErC'));
    isSide  = ~cellfun(@isempty, regexpi(roiNames, lower(side)));
    match   = hasSub & hasSes & isErC & isSide;

    if ~any(match)
        msg = sprintf('No %s ROI mask found for %s %s (filtered for ErC + %s)', ...
            upper(side), sub, ses, side);
        if failFast; error(msg); else; warning(msg); return; end
    end

    idx = find(match);
    fullPaths = fullfile(roiFolders(idx), roiNames(idx));
    [~, sortIdx] = sort(cellfun(@numel, fullPaths));
    roiPath = fullPaths{sortIdx(1)};
end


function manifest = load_run_manifest(manifestFile)
% LOAD_RUN_MANIFEST  Which run of each task the preparation actually used.
%
% Returns a containers.Map keyed 'sub|ses|task' -> 'run-2' (or 'no-run').
% An empty map when the file is absent, which is what an older GLM_runauto
% looks like; the caller then matches on subject/session/task alone.

manifest = containers.Map('KeyType', 'char', 'ValueType', 'char');

if isempty(manifestFile) || ~isfile(manifestFile)
    return;
end

fid = fopen(manifestFile, 'r');
if fid == -1
    warning('load_run_manifest:CannotRead', 'Cannot read %s', manifestFile);
    return;
end

lineNo = 0;
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    lineNo = lineNo + 1;
    if lineNo == 1, continue; end            % header
    if isempty(strtrim(line)), continue; end

    parts = strsplit(line, sprintf('\t'));
    if numel(parts) < 4, continue; end

    key = sprintf('%s|%s|%s', strtrim(parts{1}), strtrim(parts{2}), strtrim(parts{3}));
    manifest(key) = strtrim(parts{4});
end
fclose(fid);
end


function runLabel = run_entity_of(name)
% RUN_ENTITY_OF  'run-2' out of a filename fragment, or '' when it has none.
runLabel = '';
tok = regexp(name, '_run-([A-Za-z0-9]+)', 'tokens', 'once');
if ~isempty(tok)
    runLabel = ['run-' tok{1}];
end
end


function out = sort_by_trailing_number(fileList)
% Sorts a list of filenames by the number at the end.
% e.g., _00001.nii, _00002.nii, ..., _00390.nii
    n = numel(fileList);
    nums = zeros(n, 1);
    for i = 1:n
        [~, name] = fileparts(fileList{i});
        tok = regexp(name, '(\d+)$', 'tokens', 'once');
        if isempty(tok)
            nums(i) = i;
        else
            nums(i) = str2double(tok{1});
        end
    end
    [~, idx] = sort(nums);
    out = fileList(idx);
end


function write_error_log(OUTPUT_ROOT, sub, ses, errText)
% Writes an error log file for a failed session.
    try
        logDir = fullfile(OUTPUT_ROOT, '_logs');
        if ~exist(logDir, 'dir'); mkdir(logDir); end
        fid = fopen(fullfile(logDir, sprintf('%s_%s_ERROR.txt', sub, ses)), 'w');
        if fid > 0
            fprintf(fid, '%s\n', errText);
            fclose(fid);
        end
    catch
    end
end


function write_status_csv(OUTPUT_ROOT, status)
% Writes a CSV summary of which sessions succeeded or failed.
    try
        logDir = fullfile(OUTPUT_ROOT, '_logs');
        if ~exist(logDir, 'dir'); mkdir(logDir); end
        csvFile = fullfile(logDir, 'batch_status.csv');

        fid = fopen(csvFile, 'w');
        if fid < 0; return; end

        fprintf(fid, 'subject,session,ok\n');
        for i = 1:numel(status)
            fprintf(fid, '%s,%s,%d\n', status(i).sub, status(i).ses, status(i).ok);
        end
        fclose(fid);
        fprintf('Status CSV: %s\n', csvFile);
    catch
    end
end


function allowedPairs = load_allowed_pairs(listFile)
% Reads subses_list.txt and returns a containers.Map of allowed "sub-XX ses-YY" pairs.
% File format: one "sub-XX ses-YY" pair per line (space-separated).
    allowedPairs = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    fid = fopen(listFile, 'r');
    if fid < 0
        error('Could not open subject list file: %s', listFile);
    end
    cleanUp = onCleanup(@() fclose(fid));
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line) || line(1) == '#'
            continue;
        end
        % Each line is "sub-XX ses-YY"
        allowedPairs(line) = true;
    end
    fprintf('Loaded %d subject-session pairs from %s\n', allowedPairs.Count, listFile);
end


function cleanup_dir_safe(dirPath)
% Remove a directory tree if it exists; swallow any errors.
    try
        if exist(dirPath, 'dir')
            rmdir(dirPath, 's');
        end
    catch
    end
end


function gb = dir_size_gb(dirPath)
% Returns the size of a directory in GB using du.
    [~, out] = system(sprintf('du -sb "%s" 2>/dev/null', dirPath));
    tok = regexp(out, '(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        gb = str2double(tok{1}) / 1e9;
    else
        gb = NaN;
    end
end
