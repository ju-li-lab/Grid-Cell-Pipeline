function out = newest_match(folder, regexPattern)
d = dir(folder);
bestTime = -Inf;
out = '';
for i=1:numel(d)
    if d(i).isdir, continue; end
    if ~isempty(regexp(d(i).name, regexPattern, 'once'))
        if d(i).datenum > bestTime
            bestTime = d(i).datenum;
            out = fullfile(folder, d(i).name);
        end
    end
end
end

% ===================== CHOOSING WHICH RUN TO USE =====================
%
% A session can hold several runs of the same scan. The run numbers are
% counted per modality and carry no cross-modality meaning: if the subject
% climbed out of the scanner between sequences, run-2 of the T2w and run-2 of
% a task are from different visits and must not be paired.
%
% So the BOLD run being preprocessed is the anchor, and every other input is
% matched to it: an explicit choice in run_selection.tsv first, otherwise the
% run acquired closest in time, otherwise a refusal (or a loud guess).
