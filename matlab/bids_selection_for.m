function runLabel = bids_selection_for(sel, varargin)
% BIDS_SELECTION_FOR  Look one type up in a loaded run selection.
%
%   runLabel = bids_selection_for(sel, 'task-run1')
%   runLabel = bids_selection_for(sel, 'reverse', 'task-reverse', 'epi')
%
%   Types are given in order of preference; the first one present in the
%   selection wins. This is how the aliases are handled — a reverse-PE EPI may
%   be listed as 'reverse' (what scan_multirun.sh writes now), as
%   'task-reverse' (what older files say) or as 'epi' (a reverse EPI in fmap/).
%
%   Returns '' when nothing was selected, which means "decide automatically".
%
%   See also BIDS_LOAD_SELECTION, BIDS_PICK_RUN.

runLabel = '';

for i = 1:numel(varargin)
    t = varargin{i};
    if isempty(t), continue; end
    field = regexprep(t, '[^a-zA-Z0-9]', '_');
    if isfield(sel, field) && ~isempty(sel.(field))
        runLabel = strtrim(sel.(field));
        return;
    end
end
end
