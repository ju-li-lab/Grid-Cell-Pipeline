function C = bids_list_runs(folders, filt)
% BIDS_LIST_RUNS  All files in a session that match one BIDS suffix.
%
%   C = bids_list_runs(folder, filt)
%   C = bids_list_runs({folderA, folderB}, filt)
%
%   filt is a struct with any of:
%     .sub      'sub-01s13'   require this subject entity
%     .ses      'ses-01'      require this session entity
%     .suffix   'T2w'         require this BIDS suffix ('bold', 'epi', ...)
%     .task     'run1'        require this task label
%     .contains 'phasediff'   substring that must appear in the filename
%     .regexp   '...'         regular expression the filename must match
%     .ext      {'.nii','.nii.gz'}   extensions to accept (default: NIfTI)
%     .prefix   ''            required SPM prefix; '' means raw BIDS names only
%                             (omit the field to accept any prefix)
%
%   Returns a struct array, sorted by run number then by name, with fields:
%     .path    full path
%     .name    filename
%     .ent     the bids_entities struct
%     .runNum  run number, NaN when the file carries no run entity
%     .acqSec  acquisition time in seconds since midnight, NaN when unknown
%     .series  SeriesNumber from the JSON sidecar, NaN when unknown
%
%   When a session holds both foo.nii and foo.nii.gz, only the .nii is kept —
%   they are the same acquisition and SPM can only read the uncompressed one.
%
%   See also BIDS_ENTITIES, BIDS_PICK_RUN.

if ischar(folders) || isstring(folders)
    folders = {char(folders)};
end
if nargin < 2, filt = struct(); end

exts = {'.nii', '.nii.gz'};
if isfield(filt, 'ext') && ~isempty(filt.ext)
    exts = filt.ext;
    if ischar(exts), exts = {exts}; end
end

C = struct('path', {}, 'name', {}, 'ent', {}, 'runNum', {}, 'acqSec', {}, 'series', {});
seenStems = {};

for f = 1:numel(folders)
    folder = folders{f};
    if isempty(folder) || ~isfolder(folder), continue; end

    d = dir(folder);
    % Sort so .nii is seen before .nii.gz of the same stem
    [~, order] = sort({d.name});
    d = d(order);

    for i = 1:numel(d)
        if d(i).isdir, continue; end
        nm = d(i).name;
        if startsWith(nm, '._'), continue; end   % macOS resource forks

        if ~any(cellfun(@(e) endsWith(lower(nm), e), exts)), continue; end

        E = bids_entities(nm);

        if isfield(filt, 'sub')    && ~isempty(filt.sub)    && ~strcmp(E.sub, filt.sub),   continue; end
        if isfield(filt, 'ses')    && ~isempty(filt.ses)    && ~strcmp(E.ses, filt.ses),   continue; end
        if isfield(filt, 'suffix') && ~isempty(filt.suffix) && ~strcmp(E.suffix, filt.suffix), continue; end
        if isfield(filt, 'task')   && ~isempty(filt.task)   && ~strcmp(E.task, filt.task), continue; end
        if isfield(filt, 'prefix') && ~strcmp(E.prefix, filt.prefix), continue; end
        if isfield(filt, 'contains') && ~isempty(filt.contains) && ~contains(nm, filt.contains), continue; end
        if isfield(filt, 'regexp') && ~isempty(filt.regexp) && isempty(regexp(nm, filt.regexp, 'once')), continue; end

        % Same acquisition already seen as .nii — skip the .nii.gz twin
        if ismember(E.stem, seenStems), continue; end
        seenStems{end+1} = E.stem; %#ok<AGROW>

        p = fullfile(folder, nm);
        [acqSec, series] = bids_acq_time(p);

        C(end+1) = struct('path', p, 'name', nm, 'ent', E, ...
                          'runNum', E.runNum, 'acqSec', acqSec, 'series', series); %#ok<AGROW>
    end
end

if isempty(C), return; end

% Sort by run number, files without a run entity last. sort() is stable, so
% sorting by name first leaves same-run files in alphabetical order.
[~, nameOrder] = sort({C.name});
C = C(nameOrder);
runs = [C.runNum];
runs(isnan(runs)) = Inf;
[~, runOrder] = sort(runs);
C = C(runOrder);
end
