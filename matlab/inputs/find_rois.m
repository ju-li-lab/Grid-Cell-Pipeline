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
