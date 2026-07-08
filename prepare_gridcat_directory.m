function prepare_gridcat_directory(BIDS_ROOT, OUTPUT_ROOT, varargin)
% prepare_gridcat_directory - Prepare data for GridCat analysis from BIDS
%
% Usage:
%   prepare_gridcat_directory(BIDS_ROOT, OUTPUT_ROOT)
%   prepare_gridcat_directory(BIDS_ROOT, OUTPUT_ROOT, 'IncludeTasks', {'task-1','task-2','task-3'})
%   prepare_gridcat_directory(BIDS_ROOT, OUTPUT_ROOT, 'DryRun', true)
%
% Inputs:
%   BIDS_ROOT   - Path to BIDS directory
%   OUTPUT_ROOT - Path to output directory (will create GLM_runauto here)
%
% Optional Name-Value Pairs:
%   'IncludeTasks'  - Cell array of tasks to include (default: all except rest/reverse)
%   'ExcludeTasks'  - Cell array of tasks to exclude (default: {'task-rest','task-reverse'})
%   'FuncPrefix'    - Prefix for functional files (default: 'u')
%   'ROIPrefix'     - Prefix for ROI files (default: 'r')
%   'RPPrefix'      - Prefix for motion regressor files (default: 'rp')
%   'CopyMode'      - 'copy' or 'symlink' (default: 'copy')
%   'DryRun'        - If true, only shows what would be done (default: false)
%   'Verbose'       - Print detailed progress (default: true)
%
% Output Structure:
%   OUTPUT_ROOT/GLM_runauto/
%       ├── functional_scans_split/  (3D split volumes)
%       ├── rp_txt/                  (motion regressors)
%       ├── ROI/                     (bilateral masks)
%       └── _logs/                   (processing log)
%
% Example:
%   prepare_gridcat_directory('/sc-projects/.../b2_bids', '/sc-projects/.../analysis', ...
%       'IncludeTasks', {'task-1','task-2','task-3'}, 'CopyMode', 'symlink')

%% Parse inputs
p = inputParser;
addRequired(p, 'BIDS_ROOT', @(x) ischar(x) || isstring(x));
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

parse(p, BIDS_ROOT, OUTPUT_ROOT, varargin{:});
cfg = p.Results;

% Ensure paths exist
assert(isfolder(cfg.BIDS_ROOT), 'BIDS_ROOT not found: %s', cfg.BIDS_ROOT);

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

if ~cfg.DryRun
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
    fprintf(logFID, 'BIDS Root: %s\n', cfg.BIDS_ROOT);
    fprintf(logFID, 'Output Root: %s\n\n', cfg.OUTPUT_ROOT);
else
    logFID = 1; % stdout
    fprintf('\n=== DRY RUN MODE - No files will be modified ===\n\n');
end

%% Discover BIDS structure
if ~isempty(cfg.SubjectList) && isfile(cfg.SubjectList)
    fprintf('Using filtered subject list: %s\n', cfg.SubjectList);
    subjects = load_subjects_from_list(cfg.BIDS_ROOT, cfg.SubjectList);
else
    fprintf('Scanning BIDS directory...\n');
    subjects = discover_bids_subjects(cfg.BIDS_ROOT);
end

if isempty(subjects)
    error('No valid BIDS subjects found in %s', cfg.BIDS_ROOT);
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
sesItems = struct('sub', {}, 'ses', {}, 'funcDir', {}, 'anatDir', {}, 'tasks', {});
idx = 0;
for s = 1:numel(subjects)
    sub = subjects(s).name;
    for sess = 1:numel(subjects(s).sessions)
        ses = subjects(s).sessions{sess};
        sesDir = fullfile(cfg.BIDS_ROOT, sub, ses);
        fDir   = fullfile(sesDir, 'func');
        aDir   = fullfile(sesDir, 'anat');

        tasks = find_session_tasks(fDir, cfg);
        if isempty(tasks); continue; end

        idx = idx + 1;
        sesItems(idx).sub     = sub;
        sesItems(idx).ses     = ses;
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
% Each iteration is fully independent: reads from BIDS, writes to its own
% output files. No file naming conflicts because filenames include sub/ses.

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
% PROCESS_ONE_SESSION  Process a single subject-session (parfor-safe).
% Returns a struct array with one entry per task (+ ROI info on last entry).

    sub  = item.sub;
    ses  = item.ses;
    funcDir = item.funcDir;
    anatDir = item.anatDir;
    tasks   = item.tasks;
    nTasks  = numel(tasks);

    results = struct('sub', {}, 'ses', {}, 'task', {}, ...
        'func_found', {}, 'func_volumes', {}, 'rp_found', {}, ...
        'roi_left_found', {}, 'roi_right_found', {}, 'roi_bilat_created', {});

    for t = 1:nTasks
        task = tasks{t};
        r = struct('sub', sub, 'ses', ses, 'task', task, ...
            'func_found', false, 'func_volumes', 0, 'rp_found', false, ...
            'roi_left_found', false, 'roi_right_found', false, 'roi_bilat_created', false);

        % 1. Find and copy/split functional file
        funcPattern = sprintf('%s%s_%s_%s%s.nii', funcPrefix, sub, ses, task, funcSuffix);
        funcFiles = dir(fullfile(funcDir, funcPattern));

        if isempty(funcFiles)
            warning('No functional file found for %s %s %s', sub, ses, task);
        else
            funcFile = fullfile(funcDir, funcFiles(1).name);
            r.func_found = true;

            if ~dryRun
                targetFunc = fullfile(FUNC_DIR, funcFiles(1).name);
                copy_or_link(funcFile, targetFunc, copyMode, verbose);

                nVols = split_4d_to_3d(targetFunc, verbose);
                r.func_volumes = nVols;

                delete(targetFunc);  % Remove 4D copy after split
                fprintf('    %s %s %s: split into %d volumes\n', sub, ses, task, nVols);
            end
        end

        % 2. Find and copy motion regressors
        rpPattern = sprintf('%s*%s*%s*%s*.txt', rpPrefix, sub, ses, task);
        rpFiles = dir(fullfile(funcDir, rpPattern));
        rpFiles = rpFiles(~startsWith({rpFiles.name}, '._'));

        if isempty(rpFiles)
            warning('No motion regressor found for %s %s %s', sub, ses, task);
        else
            r.rp_found = true;
            if ~dryRun
                rpFile  = fullfile(funcDir, rpFiles(1).name);
                targetRP = fullfile(RP_DIR, rpFiles(1).name);
                copy_or_link(rpFile, targetRP, copyMode, verbose);
            end
        end

        results = [results, r]; %#ok<AGROW>
    end

    % 3. Process ROI masks (once per session, recorded on last task entry)
    roiAllPattern = sprintf('%s*%s*%s*.nii', roiPrefix, sub, ses);
    roiAllFiles = dir(fullfile(anatDir, roiAllPattern));
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
        targetLeft   = fullfile(ROI_DIR, roiLeftFiles(1).name);
        targetRight  = fullfile(ROI_DIR, roiRightFiles(1).name);

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


function subjects = discover_bids_subjects(bidsRoot)
% Find all valid BIDS subjects and their sessions
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
% includes subject-session pairs from the filtered list file.

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
            fprintf('  WARNING: Skipping %s %s — func/ not found\n', sub, ses);
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
% Find all tasks in a session's func directory
tasks = {};

if ~isfolder(funcDir)
    return;
end

% Find all functional files matching prefix + suffix
% FuncSuffix is '_bold' (realign_unwarp/realign_only) or '_bold_dc' (topup+applytopup)
pattern = sprintf('%s*%s.nii', cfg.FuncPrefix, cfg.FuncSuffix);
files = dir(fullfile(funcDir, pattern));
files = files(~startsWith({files.name}, '._'));

for i = 1:numel(files)
    % Extract task from filename
    tok = regexp(files(i).name, '_task-([^_]+)_', 'tokens', 'once');
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
