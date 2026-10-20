function sel = bids_load_selection(tsvFile, SUB, SES)
% BIDS_LOAD_SELECTION  Read run_selection.tsv for one subject-session.
%
%   sel = bids_load_selection('/path/run_selection.tsv', 'sub-01s13', 'ses-01')
%
%   Returns a struct whose fields are the selection types with non-identifier
%   characters replaced by '_', each holding the chosen run label:
%
%     sel.T2w          'run-1'
%     sel.task_run1    'run-2'
%     sel.task_run2    'run-2'
%     sel.reverse      'run-1'
%     sel.fieldmap     'no-run'
%
%   Column positions come from the header row when there is one, so both the
%   old five-column file (subject session type available_runs selected_run)
%   and the current one with extra columns are read correctly. Rows whose
%   selected_run is blank are skipped — a blank means "decide automatically".
%
%   Look a type up with bids_selection_for(sel, type), which also understands
%   the aliases (task-reverse / reverse, magnitude / magnitude1).
%
%   See also BIDS_PICK_RUN, BIDS_SELECTION_FOR.

sel = struct();

if isempty(tsvFile) || ~isfile(tsvFile)
    return;
end

fid = fopen(tsvFile, 'r');
if fid == -1
    warning('bids_load_selection:CannotRead', 'Cannot read run selection file: %s', tsvFile);
    return;
end
cleanupFid = onCleanup(@() fclose(fid)); %#ok<NASGU>

% Default column order = the historical five-column layout
cols = struct('subject', 1, 'session', 2, 'type', 3, 'selected_run', 5);
lineNo = 0;

while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    lineNo = lineNo + 1;

    if isempty(strtrim(line)), continue; end
    parts = strsplit(line, sprintf('\t'), 'CollapseDelimiters', false);
    parts = cellfun(@strtrim, parts, 'UniformOutput', false);

    % Header row: remember where each column actually is
    if lineNo == 1 && any(strcmpi(parts{1}, {'subject', 'sub'}))
        cols = struct();
        for i = 1:numel(parts)
            key = lower(regexprep(parts{i}, '[^a-zA-Z0-9]', '_'));
            if ~isempty(key) && isvarname(key)
                cols.(key) = i;
            end
        end
        continue;
    end

    fSub  = column_value(parts, cols, 'subject');
    fSes  = column_value(parts, cols, 'session');
    fType = column_value(parts, cols, 'type');
    fRun  = column_value(parts, cols, 'selected_run');

    if ~strcmp(fSub, SUB) || ~strcmp(fSes, SES), continue; end
    if isempty(fType) || isempty(fRun), continue; end

    field = regexprep(fType, '[^a-zA-Z0-9]', '_');
    if ~isvarname(field), continue; end
    sel.(field) = fRun;
end
end

function v = column_value(parts, cols, name)
v = '';
if ~isfield(cols, name), return; end
idx = cols.(name);
if idx >= 1 && idx <= numel(parts)
    v = parts{idx};
end
end
