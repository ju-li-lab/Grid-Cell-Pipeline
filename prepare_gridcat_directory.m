function prepare_gridcat_directory(PREPROC_ROOT, OUTPUT_ROOT, varargin)
% prepare_gridcat_directory - Flatten the preprocessing derivatives for GridCAT
%
% Usage:
%   prepare_gridcat_directory(PREPROC_ROOT, OUTPUT_ROOT)
%   prepare_gridcat_directory(PREPROC_ROOT, OUTPUT_ROOT, 'IncludeTasks', {'task-1','task-2'})
%   prepare_gridcat_directory(PREPROC_ROOT, OUTPUT_ROOT, 'DryRun', true)
%
% Inputs:
%   PREPROC_ROOT - The spm-preproc derivatives dataset written by
%                  run_spm_preproc.m (<DERIV_ROOT>/spm-preproc), laid out as
%                  sub-XX/ses-YY/{func,anat}. Not the raw BIDS directory:
%                  nothing GridCAT needs lives there any more.
%   OUTPUT_ROOT  - Path to output directory (will create GLM_runauto here)
%
% Optional Name-Value Pairs:
%   'IncludeTasks'  - Cell array of tasks to include (default: all except rest/reverse)
%   'ExcludeTasks'  - Cell array of tasks to exclude (default: {'task-rest','task-reverse'})
%   'FuncPrefix'    - Prefix for functional files (default: 'u')
%   'ROIPrefix'     - Prefix for ROI files (default: 'r')
%   'RPPrefix'      - Prefix for motion regressor files (default: 'rp')
%   'CopyMode'      - 'copy' or 'symlink' (default: 'copy')
%   'Reset'         - Empty the split/rp/ROI folders first (default: true).
%                     EventFiles/ is never touched — you put those there.
%   'DryRun'        - If true, only shows what would be done (default: false)
%   'Verbose'       - Print detailed progress (default: true)
%
% Output Structure:
%   OUTPUT_ROOT/GLM_runauto/
%       ├── functional_scans_split/  (3D split volumes)
%       ├── rp_txt/                  (motion regressors)
%       ├── ROI/                     (bilateral masks)
%       ├── EventFiles/Event_tables/ (yours — event tables, left alone)
%       ├── run_manifest.tsv         (which run of each task was used)
%       └── _logs/                   (processing log)
%
% RUN ENTITIES ARE DROPPED HERE, ON PURPOSE
% -----------------------------------------
% A session can hold several runs of a task, but exactly one of them was
% preprocessed. Inside GLM_runauto every file is therefore renamed to the
% plain <sub>_<ses>_task-<label> form, without the run entity. That is what
% makes the event tables line up: they frequently carry no run in their names,
% and matching a run-less event file against u..._task-run1_run-2_bold_0001.nii
% either fails or, worse, matches the wrong run.
%
% Nothing is lost — run_manifest.tsv records which run each flattened file came
% from, and the derivatives keep their full BIDS names.
%
% Example:
%   prepare_gridcat_directory('/sc-projects/.../derivatives/spm-preproc', ...
%       '/sc-projects/.../analysis', 'IncludeTasks', {'task-1','task-2'})

%% Parse inputs
p = inputParser;
addRequired(p, 'PREPROC_ROOT', @(x) ischar(x) || isstring(x));
addRequired(p, 'OUTPUT_ROOT', @(x) ischar(x) || isstring(x));
addParameter(p, 'IncludeTasks', {}, @iscell);
addParameter(p, 'ExcludeTasks', {'task-rest', 'task-reverse'}, @iscell);
addParameter(p, 'FuncPrefix', 'u', @ischar);
addParameter(p, 'ROIPrefix', 'r', @ischar);
addParameter(p, 'RPPrefix', 'rp', @ischar);
addParameter(p, 'CopyMode', 'copy', @(x) ismember(x, {'copy', 'symlink'}));
addParameter(p, 'DryRun', false, @islogical);
addParameter(p, 'Verbose', true, @islogical);
addParameter(p, 'SPMPath', '', @ischar);
addParameter(p, 'FuncSuffix', '_bold', @ischar);  % '_bold' or '_bold_dc' (topup)
addParameter(p, 'SubjectList', '', @ischar);       % path to subses_list.txt (filtered)
addParameter(p, 'MaxWorkers', 0, @isnumeric);      % 0 = auto (ncores-1), >0 = explicit
addParameter(p, 'Reset', true, @(x) islogical(x) || isnumeric(x));

parse(p, PREPROC_ROOT, OUTPUT_ROOT, varargin{:});
cfg = p.Results;
cfg.Reset = logical(cfg.Reset);

assert(isfolder(cfg.PREPROC_ROOT), [ ...
    'Preprocessing derivatives not found:\n  %s\n' ...
    'Run the preprocessing first (menu option 2), or check DERIV_ROOT in\n' ...
    'pipeline_config.cfg.'], cfg.PREPROC_ROOT);

% Initialize SPM if provided
if ~isempty(cfg.SPMPath)
    addpath(cfg.SPMPath);
end
try
    spm('defaults', 'fmri');
    spm_jobman('initcfg');
catch
    warning('SPM not found or not initialized. File splitting may fail.');
end

%% Create output directory structure
GLM_DIR = fullfile(cfg.OUTPUT_ROOT, 'GLM_runauto');
FUNC_DIR = fullfile(GLM_DIR, 'functional_scans_split');
RP_DIR = fullfile(GLM_DIR, 'rp_txt');
ROI_DIR = fullfile(GLM_DIR, 'ROI');
LOG_DIR = fullfile(GLM_DIR, '_logs');
MANIFEST_FILE = fullfile(GLM_DIR, 'run_manifest.tsv');

if ~cfg.DryRun
    % Start from empty folders. Preparing twice — after changing the smoothing
    % kernel, or after a different run was selected — would otherwise leave the
    % previous flattening in place, and GridCAT picks up files by prefix: the
    % stale volumes would be analysed alongside the new ones.
    if cfg.Reset
        for d = {FUNC_DIR, RP_DIR, ROI_DIR}
            if isfolder(d{1})
                fprintf('Clearing %s\n', d{1});
                rmdir(d{1}, 's');
            end
        end
    end

    if ~isfolder(cfg.OUTPUT_ROOT), mkdir(cfg.OUTPUT_ROOT); end
    if ~isfolder(GLM_DIR), mkdir(GLM_DIR); end
    if ~isfolder(FUNC_DIR), mkdir(FUNC_DIR); end
    if ~isfolder(RP_DIR), mkdir(RP_DIR); end
    if ~isfolder(ROI_DIR), mkdir(ROI_DIR); end
    if ~isfolder(LOG_DIR), mkdir(LOG_DIR); end
end

%% Initialize logging
logFile = fullfile(LOG_DIR, sprintf('preparation_log_%s.txt', datestr(now, 'yyyymmdd_HHMMSS')));
summaryFile = fullfile(LOG_DIR, 'file_summary.csv');

if ~cfg.DryRun
    logFID = fopen(logFile, 'w');
    fprintf(logFID, 'GridCat Data Preparation Log\n');
    fprintf(logFID, 'Started: %s\n', datestr(now));
    fprintf(logFID, 'Preprocessing derivatives: %s\n', cfg.PREPROC_ROOT);
    fprintf(logFID, 'Output Root: %s\n\n', cfg.OUTPUT_ROOT);
else
    logFID = 1; % stdout
    fprintf('\n=== DRY RUN MODE - No files will be modified ===\n\n');
end

%% Discover BIDS structure
if ~isempty(cfg.SubjectList) && isfile(cfg.SubjectList)
    fprintf('Using filtered subject list: %s\n', cfg.SubjectList);
    subjects = load_subjects_from_list(cfg.PREPROC_ROOT, cfg.SubjectList);
else
    fprintf('Scanning the preprocessing derivatives...\n');
    subjects = discover_bids_subjects(cfg.PREPROC_ROOT);
end

if isempty(subjects)
    error(['No preprocessed subjects found in %s\n' ...
           'Has the preprocessing run? It writes sub-XX/ses-YY/func there.'], cfg.PREPROC_ROOT);
end

fprintf('Found %d subjects\n', numel(subjects));
if logFID > 1
    fprintf(logFID, 'Found %d subjects:\n', numel(subjects));
    for i = 1:numel(subjects)
        fprintf(logFID, '  %s: %d sessions\n', subjects(i).name, numel(subjects(i).sessions));
    end
    fprintf(logFID, '\n');
end

%% Build flat list of session work items (sequential planning phase)
% Flatten subjects × sessions into a simple array for parallel dispatch.
sesItems = struct('sub', {}, 'ses', {}, 'funcDir', {}, 'anatDir', {}, 'tasks', {}, 'sesDir', {});
idx = 0;
for s = 1:numel(subjects)
    sub = subjects(s).name;
    for sess = 1:numel(subjects(s).sessions)
        ses = subjects(s).sessions{sess};
        sesDir = fullfile(cfg.PREPROC_ROOT, sub, ses);
        fDir   = fullfile(sesDir, 'func');
        aDir   = fullfile(sesDir, 'anat');

        tasks = find_session_tasks(fDir, cfg);
        if isempty(tasks); continue; end

        idx = idx + 1;
        sesItems(idx).sub     = sub;
        sesItems(idx).ses     = ses;
        sesItems(idx).sesDir  = sesDir;
        sesItems(idx).funcDir = fDir;
        sesItems(idx).anatDir = aDir;
        sesItems(idx).tasks   = tasks;
    end
end

totalSessions = numel(sesItems);
fprintf('Prepared %d session work items.\n', totalSessions);

if logFID > 1
    for i = 1:totalSessions
        fprintf(logFID, '  %s / %s: %d tasks\n', sesItems(i).sub, sesItems(i).ses, numel(sesItems(i).tasks));
    end
    fprintf(logFID, '\n');
end

%% Set up parallel pool (if available and not dry-run)
useParallel = ~cfg.DryRun && totalSessions > 1 && ...
    license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));

if useParallel
    if cfg.MaxWorkers > 0
        nWorkers = min(cfg.MaxWorkers, totalSessions);
    else
        nWorkers = min(totalSessions, feature('numcores') - 1);
        nWorkers = max(nWorkers, 1);
    end

    pool = gcp('nocreate');
    if ~isempty(pool); delete(pool); end

    fprintf('Starting parallel pool with %d workers...\n', nWorkers);
    pool = parpool('local', nWorkers);

    % Initialize SPM on all workers
    spmPath_ = char(cfg.SPMPath);
    if ~isempty(spmPath_)
        cmd = sprintf([ ...
            'addpath(''%s''); spm(''defaults'',''fmri''); spm_jobman(''initcfg''); ' ...
            'spm_get_defaults(''cmdline'', true);'], ...
            strrep(spmPath_, '''', ''''''));
        pctRunOnAll(cmd);
    end
    fprintf('Parallel pool ready.\n\n');
else
    if ~cfg.DryRun && totalSessions > 1
        fprintf('Parallel Computing Toolbox not available — processing sequentially.\n\n');
    end
end

%% Process all sessions (parallel when available)
% Each iteration is fully independent: reads one session's derivatives, writes
% to its own output files. No naming conflicts, because every output filename
% is <sub>_<ses>_task-<label> and each session is handled exactly once.

% Pre-extract scalars/strings for parfor broadcast
funcPrefix_ = cfg.FuncPrefix;
funcSuffix_ = cfg.FuncSuffix;
rpPrefix_   = cfg.RPPrefix;
roiPrefix_  = cfg.ROIPrefix;
copyMode_   = cfg.CopyMode;
verbose_    = cfg.Verbose;
dryRun_     = cfg.DryRun;

% Allocate per-session result containers
allResults = cell(totalSessions, 1);

if useParallel
    parfor k = 1:totalSessions
        allResults{k} = process_one_session(sesItems(k), ...
            FUNC_DIR, RP_DIR, ROI_DIR, ...
            funcPrefix_, funcSuffix_, rpPrefix_, roiPrefix_, ...
            copyMode_, verbose_, dryRun_);
    end
else
    for k = 1:totalSessions
        fprintf('\n[%d/%d] Processing %s / %s\n', k, totalSessions, sesItems(k).sub, sesItems(k).ses);
        allResults{k} = process_one_session(sesItems(k), ...
            FUNC_DIR, RP_DIR, ROI_DIR, ...
            funcPrefix_, funcSuffix_, rpPrefix_, roiPrefix_, ...
            copyMode_, verbose_, dryRun_);
    end
end

%% Merge results and write summary CSV
summary = [allResults{:}];

%% Write the run manifest
% One line per flattened task, saying which run of it was preprocessed. The
% GridCAT analysis reads this to pair each task with the right event table,
% and it is the only place where the run entity survives the flattening.
if ~cfg.DryRun && ~isempty(summary)
    write_run_manifest(MANIFEST_FILE, summary);
    fprintf('Run manifest written to: %s\n', MANIFEST_FILE);
end

if logFID > 1
    for i = 1:numel(summary)
        fprintf(logFID, '  %s %s %s: func=%d (%d vols), rp=%d, roi_L=%d, roi_R=%d, bilat=%d\n', ...
            summary(i).sub, summary(i).ses, summary(i).task, ...
            summary(i).func_found, summary(i).func_volumes, summary(i).rp_found, ...
            summary(i).roi_left_found, summary(i).roi_right_found, summary(i).roi_bilat_created);
    end
end

if ~cfg.DryRun && ~isempty(summary)
    write_summary_csv(summaryFile, summary);
    fprintf('\nSummary written to: %s\n', summaryFile);
end

%% Close log
if logFID > 1
    fprintf(logFID, '\nCompleted: %s\n', datestr(now));
    fclose(logFID);
    fprintf('\nLog written to: %s\n', logFile);
end

fprintf('\n=== Preparation Complete ===\n');
fprintf('Output directory: %s\n', GLM_DIR);
fprintf('Processed %d subject-sessions\n', totalSessions);

end

%% ======================== HELPER FUNCTIONS ========================

function results = process_one_session(item, FUNC_DIR, RP_DIR, ROI_DIR, ...
    funcPrefix, funcSuffix, rpPrefix, roiPrefix, copyMode, verbose, dryRun)
% PROCESS_ONE_SESSION  Flatten a single subject-session (parfor-safe).
%
% Returns a struct array with one entry per task (+ ROI info on last entry).
%
% Everything is renamed on the way out: the derivatives keep the full BIDS
% name of the run that was preprocessed (u<sub>_<ses>_task-run1_run-2_bold.nii),
% and GLM_runauto gets the run-free form (u<sub>_<ses>_task-run1_bold.nii).
% Only one run per task was preprocessed, so nothing becomes ambiguous, and the
% event tables — which often carry no run entity — line up by name.

    sub  = item.sub;
    ses  = item.ses;
    funcDir = item.funcDir;
    anatDir = item.anatDir;
    tasks   = item.tasks;
    nTasks  = numel(tasks);

    results = struct('sub', {}, 'ses', {}, 'task', {}, ...
        'func_found', {}, 'func_volumes', {}, 'rp_found', {}, ...
        'roi_left_found', {}, 'roi_right_found', {}, 'roi_bilat_created', {}, ...
        'source_run', {}, 'source_file', {});

    for t = 1:nTasks
        task = tasks{t};
        r = struct('sub', sub, 'ses', ses, 'task', task, ...
            'func_found', false, 'func_volumes', 0, 'rp_found', false, ...
            'roi_left_found', false, 'roi_right_found', false, 'roi_bilat_created', false, ...
            'source_run', 'no-run', 'source_file', '');

        % 1. The preprocessed 4D BOLD for this task.
        %    The run entity is not in the pattern: the derivatives hold exactly
        %    one run per task, whichever one the preprocessing selected.
        funcFile = find_one(funcDir, sprintf('^%s%s_%s_%s(_[^_]+)*%s\\.nii$', ...
            regexptranslate('escape', funcPrefix), regexptranslate('escape', sub), ...
            regexptranslate('escape', ses), regexptranslate('escape', task), ...
            regexptranslate('escape', funcSuffix)), sub, ses, task, 'functional file');

        if ~isempty(funcFile)
            r.func_found  = true;
            r.source_file = funcFile;
            r.source_run  = run_label_of(funcFile);

            if ~dryRun
                % Flat name, run entity dropped
                targetName = sprintf('%s%s_%s_%s%s.nii', funcPrefix, sub, ses, task, funcSuffix);
                targetFunc = fullfile(FUNC_DIR, targetName);
                copy_or_link(funcFile, targetFunc, copyMode, verbose);

                nVols = split_4d_to_3d(targetFunc, verbose);
                r.func_volumes = nVols;

                delete(targetFunc);  % Remove 4D copy after split
                fprintf('    %s %s %s (%s): split into %d volumes\n', ...
                    sub, ses, task, r.source_run, nVols);
            end
        end

        % 2. The motion regressors SPM wrote for that same file.
        %    rp_*.txt is named after the BOLD it came from, so the run entity
        %    has to be tolerated here too.
        rpFile = find_one(funcDir, sprintf('^%s_.*%s_%s_%s(_|\\.).*\\.txt$', ...
            regexptranslate('escape', strip_trailing_underscore(rpPrefix)), ...
            regexptranslate('escape', sub), regexptranslate('escape', ses), ...
            regexptranslate('escape', task)), sub, ses, task, 'motion regressor');

        if ~isempty(rpFile)
            r.rp_found = true;
            if ~dryRun
                targetRP = fullfile(RP_DIR, sprintf('%s_%s_%s_%s.txt', ...
                    strip_trailing_underscore(rpPrefix), sub, ses, task));
                copy_or_link(rpFile, targetRP, copyMode, verbose);
            end
        end

        results = [results, r]; %#ok<AGROW>
    end

    % 3. Process ROI masks (once per session, recorded on last task entry)
    roiAllFiles = dir(fullfile(anatDir, sprintf('%s*%s*%s*.nii', roiPrefix, sub, ses)));
    roiAllFiles = roiAllFiles(~startsWith({roiAllFiles.name}, '._'));

    namesAll = {roiAllFiles.name};
    isErC   = ~cellfun(@isempty, regexpi(namesAll, 'ErC'));
    isLeft  = ~cellfun(@isempty, regexpi(namesAll, 'left'));
    isRight = ~cellfun(@isempty, regexpi(namesAll, 'right'));

    roiLeftFiles  = roiAllFiles(isErC & isLeft);
    roiRightFiles = roiAllFiles(isErC & isRight);

    hasLeft  = ~isempty(roiLeftFiles);
    hasRight = ~isempty(roiRightFiles);

    if ~isempty(results)
        results(end).roi_left_found  = hasLeft;
        results(end).roi_right_found = hasRight;
    end

    if ~hasLeft || ~hasRight
        warning('Missing ROI masks for %s %s (left: %d, right: %d)', sub, ses, hasLeft, hasRight);
        if ~isempty(results)
            results(end).roi_bilat_created = false;
        end
        return;
    end

    if ~dryRun
        roiLeftFile  = fullfile(anatDir, roiLeftFiles(1).name);
        roiRightFile = fullfile(anatDir, roiRightFiles(1).name);

        % Drop the run entity here too, so the analysis finds one mask per side
        targetLeft   = fullfile(ROI_DIR, drop_run_entity(roiLeftFiles(1).name));
        targetRight  = fullfile(ROI_DIR, drop_run_entity(roiRightFiles(1).name));

        copy_or_link(roiLeftFile, targetLeft, copyMode, verbose);
        copy_or_link(roiRightFile, targetRight, copyMode, verbose);

        bilatFile = create_bilateral_mask(targetLeft, targetRight, sub, ses, ROI_DIR, roiPrefix);

        if ~isempty(results)
            results(end).roi_bilat_created = ~isempty(bilatFile);
        end
        fprintf('    %s %s: bilateral ROI created\n', sub, ses);
    else
        if ~isempty(results)
            results(end).roi_bilat_created = true;
        end
    end
end

function out = find_one(folder, pattern, sub, ses, task, what)
% FIND_ONE  The single file matching a pattern, or '' with a clear warning.
%
% More than one match means the derivatives hold two runs of the same task —
% which the preprocessing does not produce, but a directory that was never
% cleared between runs does. Picking one at random is how a session ends up
% analysed with the wrong run, so this refuses to choose.

out = '';
d = dir(folder);
hits = {};
for i = 1:numel(d)
    if d(i).isdir, continue; end
    if startsWith(d(i).name, '._'), continue; end
    if ~isempty(regexp(d(i).name, pattern, 'once'))
        hits{end+1} = d(i).name; %#ok<AGROW>
    end
end

if isempty(hits)
    warning('prepare_gridcat:MissingFile', 'No %s found for %s %s %s in %s', ...
            what, sub, ses, task, folder);
    return;
end

if numel(hits) > 1
    warning('prepare_gridcat:AmbiguousFile', [ ...
        'Several %s files for %s %s %s in %s:\n  %s\n' ...
        'The preprocessing writes one run per task, so this directory holds\n' ...
        'output from more than one run. Re-run the preprocessing with\n' ...
        'DERIV_RESET=auto (the default) to clear it, and skip this session.'], ...
        what, sub, ses, task, folder, strjoin(hits, sprintf('\n  ')));
    return;
end

out = fullfile(folder, hits{1});
end

function lbl = run_label_of(f)
% The run entity of a filename, or 'no-run'.
[~, base] = fileparts(f);
tok = regexp(base, '_run-([A-Za-z0-9]+)', 'tokens', 'once');
if isempty(tok)
    lbl = 'no-run';
else
    lbl = ['run-' tok{1}];
end
end

function out = drop_run_entity(name)
% Remove the _run-N entity from a BIDS filename.
out = regexprep(name, '_run-[A-Za-z0-9]+', '');
end

function p = strip_trailing_underscore(p)
p = regexprep(p, '_+$', '');
end

function write_run_manifest(manifestFile, summary)
% WRITE_RUN_MANIFEST  Which run of each task ended up in GLM_runauto.
%
% run_gridcat_analysis.m reads this to decide whether an event table named
% ..._run-2_EventData.txt belongs to the data that was actually preprocessed.
% An event file with no run entity matches whatever is listed here.

fid = fopen(manifestFile, 'w');
if fid < 0
    warning('prepare_gridcat:ManifestFailed', 'Could not write %s', manifestFile);
    return;
end

fprintf(fid, 'subject\tsession\ttask\tselected_run\tsource_file\n');
for i = 1:numel(summary)
    if ~summary(i).func_found, continue; end
    fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', summary(i).sub, summary(i).ses, ...
        summary(i).task, summary(i).source_run, summary(i).source_file);
end
fclose(fid);
end


function subjects = discover_bids_subjects(bidsRoot)
% Find all subjects and sessions in a BIDS-style tree (raw or derivative)
subjects = struct('name', {}, 'sessions', {});

d = dir(bidsRoot);
subDirs = d([d.isdir] & startsWith({d.name}, 'sub-'));

for i = 1:numel(subDirs)
    subName = subDirs(i).name;
    subPath = fullfile(bidsRoot, subName);
    
    % Find sessions
    sd = dir(subPath);
    sesDirs = sd([sd.isdir] & startsWith({sd.name}, 'ses-'));
    
    if isempty(sesDirs)
        % No session folders - might be single-session study
        % Check if func folder exists directly under subject
        if isfolder(fullfile(subPath, 'func'))
            sessions = {''}; % Empty string means no session subfolder
        else
            continue;
        end
    else
        sessions = {sesDirs.name};
    end
    
    subjects(end+1).name = subName; %#ok<AGROW>
    subjects(end).sessions = sessions;
end
end

function subjects = load_subjects_from_list(bidsRoot, listFile)
% LOAD_SUBJECTS_FROM_LIST  Read subses_list.txt and build subjects struct.
%
% Returns the same struct format as discover_bids_subjects, but only
% includes subject-session pairs from the filtered list file. bidsRoot is the
% preprocessing derivatives here, so a pair whose func/ is missing simply has
% not been preprocessed yet.

subjects = struct('name', {}, 'sessions', {});

fid = fopen(listFile, 'r');
if fid == -1
    warning('Cannot read subject list: %s', listFile);
    return;
end

% Read all lines
subSesMap = containers.Map();
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    line = strtrim(line);
    if isempty(line), continue; end

    parts = strsplit(line);
    sub = parts{1};
    if numel(parts) >= 2
        ses = parts{2};
    else
        ses = '';
    end

    % Group sessions by subject
    if subSesMap.isKey(sub)
        subSesMap(sub) = [subSesMap(sub), {ses}];
    else
        subSesMap(sub) = {ses};
    end
end
fclose(fid);

% Build subjects struct
keys = subSesMap.keys();
for i = 1:numel(keys)
    sub = keys{i};
    sessions = subSesMap(sub);
    % Verify directories exist
    validSessions = {};
    for j = 1:numel(sessions)
        ses = sessions{j};
        if isempty(ses)
            checkDir = fullfile(bidsRoot, sub, 'func');
        else
            checkDir = fullfile(bidsRoot, sub, ses, 'func');
        end
        if isfolder(checkDir)
            validSessions{end+1} = ses; %#ok<AGROW>
        else
            fprintf('  WARNING: Skipping %s %s — no preprocessed func/ in %s\n', ...
                    sub, ses, bidsRoot);
        end
    end
    if ~isempty(validSessions)
        subjects(end+1).name = sub; %#ok<AGROW>
        subjects(end).sessions = validSessions;
    end
end

fprintf('  Loaded %d subjects (%d sessions) from %s\n', ...
    numel(subjects), sum(cellfun(@numel, {subjects.sessions})), listFile);
end

function tasks = find_session_tasks(funcDir, cfg)
% FIND_SESSION_TASKS  The task labels this session has preprocessed data for.
tasks = {};

if ~isfolder(funcDir)
    return;
end

% All functional files matching prefix + suffix.
% FuncSuffix is '_bold' (realign_unwarp/realign_only) or '_bold_dc' (topup+applytopup).
% The 'sub-' after the prefix matters: with FuncPrefix='su', a bare 'su*_bold.nii'
% also matches the raw sub-XX_..._bold.nii the preprocessing staged, because a
% BIDS name starts with 'sub-'.
pattern = sprintf('%ssub-*%s.nii', cfg.FuncPrefix, cfg.FuncSuffix);
files = dir(fullfile(funcDir, pattern));
files = files(~startsWith({files.name}, '._'));

for i = 1:numel(files)
    name = files(i).name;

    % The per-task and session mean images match the same prefix pattern
    % (meanu<sub>_..._task-run1_bold.nii). They are not runs and must not
    % become tasks — this is why 'mean' had to be listed in EXCLUDE_TASKS.
    if ~isempty(regexp(name, '^mean', 'once')) || ...
       ~isempty(regexp(name, ['^' regexptranslate('escape', cfg.FuncPrefix) 'mean'], 'once'))
        continue;
    end

    % Extract the task label. The run entity, when there is one, sits between
    % the task and the suffix and is not part of the label.
    tok = regexp(name, '_task-([A-Za-z0-9]+)', 'tokens', 'once');
    if isempty(tok)
        continue;
    end

    task = ['task-' tok{1}];

    % Apply task filtering
    if ~isempty(cfg.IncludeTasks)
        if ~any(strcmp(task, cfg.IncludeTasks))
            continue;
        end
    end

    if any(strcmp(task, cfg.ExcludeTasks))
        continue;
    end

    % Add to list if not already there
    if ~any(strcmp(task, tasks))
        tasks{end+1} = task; %#ok<AGROW>
    end
end

tasks = sort(tasks);
end

function copy_or_link(source, target, mode, verbose)
% Copy or create symbolic link
if strcmp(mode, 'copy')
    copyfile(source, target);
    if verbose
        fprintf('    Copied: %s\n', source);
    end
else
    % Create symbolic link (Unix/Linux/Mac)
    [status, msg] = system(sprintf('ln -s "%s" "%s"', source, target));
    if status ~= 0
        warning('Failed to create symlink: %s. Falling back to copy.', msg);
        copyfile(source, target);
    elseif verbose
        fprintf('    Linked: %s\n', source);
    end
end
end

function nVols = split_4d_to_3d(fourDFile, verbose)
% Split 4D NIfTI into 3D volumes using SPM
try
    % Get file info
    V = spm_vol(fourDFile);
    nVols = numel(V);
    
    if nVols == 1
        % Already 3D
        if verbose
            fprintf('    File is already 3D, no splitting needed\n');
        end
        return;
    end
    
    % Split using SPM
    spm_file_split(fourDFile);
    
catch ME
    warning('Failed to split 4D file: %s\nError: %s', fourDFile, ME.message);
    nVols = 0;
end
end

function bilatFile = create_bilateral_mask(leftFile, rightFile, sub, ses, outputDir, prefix)
% Create bilateral mask by combining left and right ROIs
bilatFile = '';

try
    % Output filename
    [~, leftName] = fileparts(leftFile);
    bilatName = strrep(leftName, 'left', 'bilat');
    bilatFile = fullfile(outputDir, [bilatName '.nii']);
    
    % Use SPM imcalc to combine masks
    matlabbatch = {};
    matlabbatch{1}.spm.util.imcalc.input = {
        leftFile
        rightFile
    };
    matlabbatch{1}.spm.util.imcalc.output = bilatName;
    matlabbatch{1}.spm.util.imcalc.outdir = {outputDir};
    matlabbatch{1}.spm.util.imcalc.expression = '(i1 + i2) > 0'; % Union of masks
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 0; % Nearest neighbor
    matlabbatch{1}.spm.util.imcalc.options.dtype = 2; % uint8
    
    spm_jobman('run', matlabbatch);
    
catch ME
    warning('Failed to create bilateral mask: %s', ME.message);
    bilatFile = '';
end
end

function write_summary_csv(filename, summary)
% Write summary to CSV file
fid = fopen(filename, 'w');
if fid < 0
    warning('Could not write summary CSV');
    return;
end

% Header
fprintf(fid, 'subject,session,task,func_found,func_volumes,rp_found,roi_left,roi_right,roi_bilat\n');

% Data
for i = 1:numel(summary)
    fprintf(fid, '%s,%s,%s,%d,%d,%d,%d,%d,%d\n', ...
        summary(i).sub, summary(i).ses, summary(i).task, ...
        summary(i).func_found, summary(i).func_volumes, summary(i).rp_found, ...
        summary(i).roi_left_found, summary(i).roi_right_found, summary(i).roi_bilat_created);
end

fclose(fid);
end
