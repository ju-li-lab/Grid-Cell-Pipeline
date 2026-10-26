function outFiles = slice_time_correct(boldFiles, CFG)
% SLICE_TIME_CORRECT  Shift every slice to a common point in the TR.
%
%   outFiles = slice_time_correct(boldFiles, CFG)
%
% Returns the corrected files (the same names with SLICETIME_PREFIX in
% front). Runs with no BOLD are passed through as empty.
%
% WHY THE TIMES AND NOT THE ORDER
% -------------------------------
% SPM accepts either a slice ORDER (1:30, 30:-1:1, an interleaved vector)
% or the actual slice TIMES. This uses the times, read from each run's own
% JSON sidecar, because an order has to be guessed and a guess here is
% expensive: correcting with the wrong order is worse than not correcting
% at all. The times come from the scanner.
%
% For this dataset the protocol sheet says "Multi-Slice Mode: Interleaved"
% while the DICOM says sequential descending (30 -> 1). The DICOM wins -
% it records what happened rather than what was programmed - and reading
% the JSON means the question never has to be answered by hand.
%
% SPM switches into slice-times mode on its own when any value in `so`
% exceeds nslices, and then `refslice` is a TIME in the same units (ms),
% not a slice index. `ta` is ignored and must be 0.
%
% THE REFERENCE POINT MATTERS DOWNSTREAM
% --------------------------------------
% Everything is shifted to CFG.spm.slicetime.ref_fraction of the TR.
% At 0.5 (the middle) the data then matches SPM's microtime default of
% T0 = T/2, which is what GridCAT's microtimeOnset/microtimeResolution
% express. Change one and you must change the other, so this function
% checks the pair and says so if they disagree.

nTasks   = numel(boldFiles);
outFiles = cell(nTasks, 1);

frac = 0.5;
if isfield(CFG, 'spm') && isfield(CFG.spm, 'slicetime') && ...
   isfield(CFG.spm.slicetime, 'ref_fraction')
    frac = CFG.spm.slicetime.ref_fraction;
end
assert(frac >= 0 && frac <= 1, ...
       'slice_time_correct: ref_fraction must be between 0 and 1 (got %g).', frac);

check_microtime(CFG, frac);

for ii = 1:nTasks
    if isempty(boldFiles{ii})
        continue;
    end

    [st_s, nslices] = read_slice_timing(boldFiles{ii});
    tr = CFG.TR;

    % Slices must fit inside one TR. A SliceTiming that runs past it means
    % the JSON and the TR disagree, and the correction would be nonsense.
    assert(max(st_s) < tr, [ ...
        'slice_time_correct: SliceTiming runs to %.4f s but TR is %.4f s in %s.\n' ...
        'The sidecar and the configured TR disagree - fix one before correcting.'], ...
        max(st_s), tr, boldFiles{ii});

    refTime_ms = frac * tr * 1000;

    matlabbatch = {};
    matlabbatch{1}.spm.temporal.st.scans    = {expand_4d(boldFiles{ii})};
    matlabbatch{1}.spm.temporal.st.nslices  = nslices;
    matlabbatch{1}.spm.temporal.st.tr       = tr;
    matlabbatch{1}.spm.temporal.st.ta       = 0;            % 0 when times are given
    matlabbatch{1}.spm.temporal.st.so       = st_s(:)' * 1000;   % ms
    matlabbatch{1}.spm.temporal.st.refslice = refTime_ms;        % ms, not an index
    matlabbatch{1}.spm.temporal.st.prefix   = CFG.SLICETIME_PREFIX;

    fprintf('  %s\n', shortname(boldFiles{ii}));
    fprintf('    %d slices, acquired %s, reference %.0f ms (%.2f x TR)\n', ...
            nslices, describe_order(st_s), refTime_ms, frac);

    spm_jobman('run', matlabbatch);

    outFiles{ii} = prefixed(boldFiles{ii}, CFG.SLICETIME_PREFIX);
    assert(isfile(outFiles{ii}), ...
           'slice_time_correct: SPM did not write the expected file: %s', outFiles{ii});
end
end


% ---------------------------------------------------------------------------

function [st_s, nslices] = read_slice_timing(boldFile)
% SliceTiming in seconds, from the sidecar that bids_stage copied along.
jsonFile = regexprep(boldFile, '\.nii(\.gz)?$', '.json');
assert(isfile(jsonFile), [ ...
    'slice_time_correct: no JSON sidecar for\n  %s\n' ...
    'Slice timing needs SliceTiming from the sidecar. Either the staging\n' ...
    'did not copy it, or the BIDS conversion did not write it.'], boldFile);

j = read_json_file(jsonFile);
assert(isfield(j, 'SliceTiming') && ~isempty(j.SliceTiming), [ ...
    'slice_time_correct: %s has no SliceTiming.\n' ...
    'Without it the slice order would have to be guessed, and guessing\n' ...
    'wrong makes the data worse. Set DO_SLICE_TIMING=0, or fix the\n' ...
    'conversion so the sidecars carry SliceTiming.'], jsonFile);

st_s    = double(j.SliceTiming(:));
nslices = numel(st_s);
end


function s = describe_order(st_s)
% Name the pattern, so the log says what was assumed rather than hiding it.
[~, order] = sort(st_s);
n = numel(st_s);
if isequal(order(:)', 1:n)
    s = 'ascending (1 -> end)';
elseif isequal(order(:)', n:-1:1)
    s = 'descending (end -> 1)';
elseif all(mod(order(1:ceil(n/2)), 2) == mod(order(1), 2))
    s = sprintf('interleaved (%s first)', tern(mod(order(1),2)==0, 'even', 'odd'));
else
    s = 'irregular - taken straight from the sidecar';
end
end


function check_microtime(CFG, frac)
% SPM's microtime onset says where in the TR the design matrix is sampled.
% After correcting to `frac` of the TR, T0 should be frac*T. GridCAT holds
% these as microtimeOnset/microtimeResolution, so a mismatch here quietly
% costs sensitivity in the GLM.
if ~isfield(CFG, 'spm') || ~isfield(CFG.spm, 'GLM'), return; end
G = CFG.spm.GLM;
if ~isfield(G, 'microtimeOnset') || ~isfield(G, 'microtimeResolution'), return; end
want = round(frac * G.microtimeResolution);
if G.microtimeOnset ~= want
    warning('slice_time_correct:MicrotimeMismatch', [ ...
        'Slice timing references %.2f x TR, so microtimeOnset should be\n' ...
        '%d (= %.2f x microtimeResolution %d), but spm_options.m has %d.\n' ...
        'The GLM will model the response at the wrong point in the TR.'], ...
        frac, want, frac, G.microtimeResolution, G.microtimeOnset);
else
    fprintf('    microtimeOnset %d matches the reference point.\n', G.microtimeOnset);
end
end


function p = prefixed(f, prefix)
[d, n, e] = fileparts(f);
p = fullfile(d, [prefix n e]);
end


function s = shortname(f)
[~, n, e] = fileparts(f);
s = [n e];
end


function out = tern(cond, a, b)
if cond, out = a; else, out = b; end
end
