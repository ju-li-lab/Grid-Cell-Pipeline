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
