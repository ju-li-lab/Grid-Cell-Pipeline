function [picks, anchor] = select_task_bolds(rawFuncDir, SUB, SES, CFG, runSel, required)
% SELECT_TASK_BOLDS  One BOLD run per task label, and the anchor time.
%
% Returns a struct array with one entry per CFG.TASK_LABELS:
%   .label .path .runLabel .acqSec .reason
% and the anchor, which is what the anatomy and fieldmaps are matched against:
%   .secs  every chosen functional run
%   .sec   the earliest of them, and .name the task it belongs to
%
% Tasks are resolved in two passes. A task with only one run needs no decision
% and pins down where the functional block is; those go first, and the tasks
% that do need a decision are then matched against them. Resolving in config
% order instead would leave the first task with nothing to match to — exactly
% the case where a session has two runs of task-run1 and one of task-run2.

nTasks = numel(CFG.TASK_LABELS);
picks = repmat(struct('label', '', 'path', '', 'runLabel', '', ...
                      'acqSec', NaN, 'reason', ''), nTasks, 1);

anchor = struct('sec', NaN, 'secs', NaN, 'name', 'the functional data');

% --------- Gather the candidates for every task ----------
cands = cell(nTasks, 1);
for ii = 1:nTasks
    tLabel = CFG.TASK_LABELS{ii};
    picks(ii).label = tLabel;

    cands{ii} = bids_list_runs(rawFuncDir, struct('sub', SUB, 'ses', SES, ...
                                                  'suffix', 'bold', 'task', tLabel, ...
                                                  'prefix', ''));

    if isempty(cands{ii}) && required
        error('run_spm_preproc:NoBold', [ ...
            'No BOLD file for task-%s in %s\n' ...
            'Expected something like %s_%s_task-%s[_run-N]_bold.nii[.gz]\n' ...
            'Check the TASKS setting in config.sh against the filenames.'], ...
            tLabel, rawFuncDir, SUB, SES, tLabel);
    end
end

% --------- Pass 1: the tasks that decide themselves ----------
% Either only one run exists, or run_selection.tsv already names one.
resolved = false(nTasks, 1);
anchorSecs = [];

for ii = 1:nTasks
    if isempty(cands{ii}), continue; end
    tLabel = picks(ii).label;
    selected = bids_selection_for(runSel, ['task-' tLabel], tLabel);

    if numel(cands{ii}) > 1 && isempty(selected)
        continue;   % needs an anchor; pass 2
    end

    [p, info] = bids_pick_run(cands{ii}, struct( ...
        'label',    sprintf('BOLD task-%s', tLabel), ...
        'selected', selected, ...
        'gapSec',   CFG.RUN_MATCH_GAP_SEC, ...
        'strict',   CFG.STRICT_RUN_MATCHING));

    picks(ii) = store_pick(picks(ii), p, info);
    resolved(ii) = true;
    if ~isnan(info.acqSec), anchorSecs(end+1) = info.acqSec; end %#ok<AGROW>
    fprintf('  BOLD task-%-8s %-10s %s\n', tLabel, info.runLabel, info.reason);
end

% --------- Pass 2: the rest, matched to what pass 1 settled ----------
for ii = 1:nTasks
    if resolved(ii) || isempty(cands{ii}), continue; end
    tLabel = picks(ii).label;

    [p, info] = bids_pick_run(cands{ii}, struct( ...
        'label',      sprintf('BOLD task-%s', tLabel), ...
        'selected',   '', ...
        'anchorSec',  anchor_vector(anchorSecs), ...
        'anchorName', 'the other functional runs', ...
        'gapSec',     CFG.RUN_MATCH_GAP_SEC, ...
        'strict',     CFG.STRICT_RUN_MATCHING));

    picks(ii) = store_pick(picks(ii), p, info);
    if ~isnan(info.acqSec), anchorSecs(end+1) = info.acqSec; end %#ok<AGROW>
    fprintf('  BOLD task-%-8s %-10s %s\n', tLabel, info.runLabel, info.reason);
end

% --------- The anchor everything else is matched against ----------
secs = [picks.acqSec];
secs = secs(~isnan(secs));
if ~isempty(secs)
    anchor.secs = secs;
    anchor.sec  = min(secs);
    idx = 1;
    for ii = 1:nTasks
        if ~isnan(picks(ii).acqSec) && picks(ii).acqSec == anchor.sec
            idx = ii;
            break;
        end
    end
    anchor.name = sprintf('task-%s %s', picks(idx).label, picks(idx).runLabel);

    % Functional runs spread over more than one visit: the anatomy can only
    % match one of them, so say so before anything is coregistered.
    if max(secs) - min(secs) > CFG.RUN_MATCH_GAP_SEC
        warning('run_spm_preproc:TasksSpanBlocks', [ ...
            'The chosen functional runs span %.0f minutes, more than\n' ...
            'RUN_MATCH_GAP_MIN. They were probably not acquired in one visit to the\n' ...
            'scanner, so one session mean covers two head positions and the T2w can\n' ...
            'only match one of them. Check run_selection.tsv for %s / %s.'], ...
            (max(secs) - min(secs))/60, SUB, SES);
    end
end
end

function pick = store_pick(pick, p, info)
pick.path     = p;
pick.runLabel = info.runLabel;
pick.acqSec   = info.acqSec;
pick.reason   = info.reason;
end

function v = anchor_vector(secs)
if isempty(secs)
    v = NaN;
else
    v = secs;
end
end
