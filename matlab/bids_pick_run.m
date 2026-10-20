function [pick, info] = bids_pick_run(cands, opts)
% BIDS_PICK_RUN  Choose one run out of several, and say why.
%
%   [pick, info] = bids_pick_run(cands, opts)
%
%   cands  struct array from bids_list_runs
%   opts   struct with:
%     .label      what this is, for messages ('T2w', 'BOLD task-run1', ...)
%     .selected   run label from run_selection.tsv ('run-2', 'no-run', '')
%     .anchorSec  acquisition time(s) to match against, NaN for none. A vector
%                 is allowed and is the normal case: the distance that matters
%                 is to the NEAREST functional run, not to the first one. In a
%                 long protocol the last scan is naturally far from the first
%                 while still belonging to the same visit.
%     .anchorName what the anchor is, for messages
%     .gapSec     how far apart two scans may be and still count as the same
%                 block (default 20 min)
%     .strict     true = refuse to guess when the choice is ambiguous
%     .required   true = error when there is no candidate at all
%
%   info  .path, .runLabel, .acqSec, .gapSec, .reason, .matched
%
%   THE PROBLEM THIS SOLVES
%   -----------------------
%   A session can hold several runs of the same scan because the subject got
%   out of the scanner and came back. The run numbers are then per modality
%   and say nothing about which scans belong together: run-2 of the T2w may
%   have been acquired in the first block and run-2 of a task in the second.
%   Pairing by run number quietly coregisters the wrong anatomy.
%
%   So the order of preference is:
%     1. An explicit choice in run_selection.tsv — the human decided.
%     2. The candidate acquired closest in time to the anchor (the BOLD run
%        being preprocessed), as long as it is within .gapSec of it. This is
%        the "same visit to the scanner" test.
%     3. Nothing else is safe: with .strict the pick fails and asks for a
%        selection, otherwise it falls back to the highest run number and says
%        loudly that it guessed.
%
%   See also BIDS_LIST_RUNS, BIDS_LOAD_SELECTION.

if nargin < 2, opts = struct(); end
opts = fill_defaults(opts);

info = struct('path', '', 'runLabel', '', 'acqSec', NaN, 'gapSec', NaN, ...
              'reason', '', 'matched', true);
pick = '';

% ---- Nothing found ----------------------------------------------------------
if isempty(cands)
    if opts.required
        error('bids_pick_run:NoCandidate', 'No %s found.', opts.label);
    end
    info.reason = 'none found';
    info.matched = false;
    return;
end

% ---- Exactly one: no decision to make ---------------------------------------
if numel(cands) == 1
    pick = cands(1).path;
    info = fill_info(info, cands(1), 'only one available');
    % An explicit choice that points elsewhere is still worth flagging
    if ~isempty(opts.selected) && ~strcmp(info.runLabel, normalise_run(opts.selected))
        warning('bids_pick_run:SelectionIgnored', ...
            ['run_selection.tsv asks for %s for %s, but the only file present is %s.\n' ...
             '  Using: %s'], opts.selected, opts.label, info.runLabel, cands(1).name);
    end
    return;
end

% ---- An explicit selection wins ---------------------------------------------
if ~isempty(opts.selected)
    want = normalise_run(opts.selected);
    for i = 1:numel(cands)
        if strcmp(run_label(cands(i)), want)
            pick = cands(i).path;
            info = fill_info(info, cands(i), 'chosen in run_selection.tsv');
            info.gapSec = gap_to_anchor(cands(i), opts);
            warn_if_far(info, opts);
            return;
        end
    end
    msg = sprintf(['run_selection.tsv asks for %s for %s, but no such file exists.\n' ...
                   '  Available: %s'], opts.selected, opts.label, available_str(cands));
    if opts.strict
        error('bids_pick_run:SelectionNotFound', '%s', msg);
    end
    warning('bids_pick_run:SelectionNotFound', '%s\n  Falling back to matching by acquisition time.', msg);
end

% ---- Match by acquisition time to the anchor --------------------------------
if any(~isnan(opts.anchorSec))
    gaps = arrayfun(@(c) nearest_gap(c.acqSec, opts.anchorSec), cands);
    known = ~isnan(gaps);
    if any(known)
        candidateGaps = gaps;
        candidateGaps(~known) = Inf;
        [bestGap, idx] = min(candidateGaps);

        if bestGap <= opts.gapSec
            pick = cands(idx).path;
            info = fill_info(info, cands(idx), sprintf('acquired %s from %s', ...
                format_gap(bestGap), opts.anchorName));
            info.gapSec = bestGap;
            return;
        end

        % Every candidate is in a different block than the anchor. Saying so is
        % more useful than silently taking the nearest one.
        msg = sprintf([ ...
            'None of the %s runs was acquired near %s (closest is %s away).\n' ...
            '  This session looks like the subject left the scanner between sequences.\n' ...
            '  Available: %s\n' ...
            '  Pick one explicitly in run_selection.tsv (bash scan_multirun.sh writes it).'], ...
            opts.label, opts.anchorName, format_gap(bestGap), available_str(cands));
        if opts.strict
            error('bids_pick_run:BlockMismatch', '%s', msg);
        end
        warning('bids_pick_run:BlockMismatch', '%s\n  Using the closest one anyway.', msg);
        pick = cands(idx).path;
        info = fill_info(info, cands(idx), sprintf('closest in time (%s away — NOT the same block)', format_gap(bestGap)));
        info.gapSec = bestGap;
        info.matched = false;
        return;
    end
end

% ---- Nothing to go on -------------------------------------------------------
msg = sprintf([ ...
    'Several %s runs exist and nothing says which one to use.\n' ...
    '  Available: %s\n' ...
    '  No acquisition time in the JSON sidecars, so the runs cannot be matched\n' ...
    '  to the functional data automatically.\n' ...
    '  Run  bash scan_multirun.sh  and fill in run_selection.tsv.'], ...
    opts.label, available_str(cands));

if opts.strict
    error('bids_pick_run:Ambiguous', '%s', msg);
end

idx = numel(cands);   % bids_list_runs sorts by run number, so this is the last
pick = cands(idx).path;
info = fill_info(info, cands(idx), 'guessed: highest run number');
info.matched = false;
warning('bids_pick_run:Ambiguous', '%s\n  Guessing the highest run number: %s', msg, cands(idx).name);
end

% ============================== helpers ==============================

function opts = fill_defaults(opts)
% Only absent fields get a default — an explicitly passed false or '' stands.
defaults = struct('label', 'file', 'selected', '', 'anchorSec', NaN, ...
                  'anchorName', 'the functional data', 'gapSec', 20*60, ...
                  'strict', true, 'required', true);
fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i})
        opts.(fn{i}) = defaults.(fn{i});
    end
end
opts.strict   = logical(opts.strict);
opts.required = logical(opts.required);
end

function lbl = run_label(c)
if isempty(c.ent.run)
    lbl = 'no-run';
else
    lbl = c.ent.run;
end
end

function want = normalise_run(sel)
sel = strtrim(sel);
if any(strcmpi(sel, {'no-run', 'norun', 'none', '-'}))
    want = 'no-run';
elseif ~isempty(regexp(sel, '^\d+$', 'once'))
    want = ['run-' sel];
else
    want = sel;
end
end

function info = fill_info(info, c, reason)
info.path     = c.path;
info.runLabel = run_label(c);
info.acqSec   = c.acqSec;
info.reason   = reason;
end

function g = gap_to_anchor(c, opts)
g = nearest_gap(c.acqSec, opts.anchorSec);
end

function g = nearest_gap(candSec, anchorSecs)
% NEAREST_GAP  Distance from a candidate to the closest anchor time.
g = NaN;
if isnan(candSec), return; end
anchorSecs = anchorSecs(~isnan(anchorSecs));
if isempty(anchorSecs), return; end
g = min(abs(anchorSecs - candSec));
end

function warn_if_far(info, opts)
if ~isnan(info.gapSec) && info.gapSec > opts.gapSec
    warning('bids_pick_run:SelectedFarFromAnchor', ...
        ['The %s chosen in run_selection.tsv (%s) was acquired %s from %s.\n' ...
         '  That is more than RUN_MATCH_GAP_MIN, so they may be from different\n' ...
         '  visits to the scanner. Using it because you asked for it.'], ...
        opts.label, info.runLabel, format_gap(info.gapSec), opts.anchorName);
end
end

function s = available_str(cands)
labels = cell(1, numel(cands));
for i = 1:numel(cands)
    stamp = '';
    if ~isnan(cands(i).acqSec)
        stamp = sprintf(' @%s', clock_str(cands(i).acqSec));
    end
    labels{i} = sprintf('%s%s', run_label(cands(i)), stamp);
end
s = strjoin(labels, ', ');
end

function s = clock_str(sec)
s = sprintf('%02d:%02d:%02d', floor(sec/3600), floor(mod(sec,3600)/60), floor(mod(sec,60)));
end

function s = format_gap(sec)
if isnan(sec)
    s = 'an unknown time';
elseif sec < 90
    s = sprintf('%.0f s', sec);
elseif sec < 5400
    s = sprintf('%.0f min', sec/60);
else
    s = sprintf('%.1f h', sec/3600);
end
end
