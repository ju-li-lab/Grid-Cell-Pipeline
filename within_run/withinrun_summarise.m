function S = withinrun_summarise(csvFiles, varargin)
% WITHINRUN_SUMMARISE  Turn per-event tables into within-run views.
%
% Reads one or more CSVs written by withinrun_extract_events and answers the
% question that motivated this work: does the aligned-vs-misaligned response
% change over the course of a run? It bins events by their position within
% the run and computes, per bin, the mean aligned response, the mean
% misaligned response, and their difference (the within-run analogue of the
% single GridCAT contrast value).
%
% USAGE
%   S = withinrun_summarise('sub-01s13_ses-01_bilat_events.csv');
%   S = withinrun_summarise({fileA, fileB, ...}, 'nBins', 4, 'useSet','glm2');
%
% OPTIONS
%   'useSet'  'glm2' (default) | 'glm1' | 'all'
%             which events to include (glm2 = the test set you contrast).
%   'nBins'   number of within-run time bins (default 2 = early/late).
%             Fewer bins = more events per estimate = more power. Start at 2.
%   'plot'    true (default) -> draw figures; false -> just return S.
%
% RETURNS struct S with, per bin:
%   S.binCentre        within-run fraction at the bin centre
%   S.aligned_mean     mean response of aligned events in that bin
%   S.misaligned_mean  mean response of misaligned events
%   S.contrast         aligned_mean - misaligned_mean  (the key quantity)
%   S.contrast_sem     between-run SEM of the per-run contrast in that bin
%   S.nAligned/nMisaligned  event counts (watch these on small data!)
%
% This is intentionally a few lines of arithmetic on the event table. Once
% you trust the picture, the same table supports a proper mixed-effects
% model (response ~ alignment * within_run_frac + (1|run)) without re-running
% any GLM.

p = inputParser;
p.addParameter('useSet', 'glm2');
p.addParameter('nBins', 2);
p.addParameter('plot', true);
p.parse(varargin{:});
opt = p.Results;

if ischar(csvFiles) || isstring(csvFiles); csvFiles = cellstr(csvFiles); end

% ---- load and concatenate all event tables
T = table();
for i = 1:numel(csvFiles)
    Ti = readtable(csvFiles{i});
    T  = [T; Ti]; %#ok<AGROW>
end

% ---- restrict to the requested event set
switch lower(opt.useSet)
    case 'glm2', T = T(T.in_glm2 == 1, :);
    case 'glm1', T = T(T.in_glm1 == 1, :);
    case 'all'   % keep everything
    otherwise, error('useSet must be glm2 | glm1 | all');
end
assert(~isempty(T), 'No events left after filtering (useSet=%s).', opt.useSet);

% ---- assign each event to a within-run time bin
edges = linspace(0, 1, opt.nBins + 1);
edges(end) = edges(end) + eps;                 % include frac == 1
bin = discretize(T.within_run_frac, edges);

% unique run identity (sub|ses|run) so we can get a between-run SEM
runKey = strcat(string(T.sub),'|',string(T.ses),'|',string(T.run));
uRuns  = unique(runKey);

S.binCentre       = (edges(1:end-1) + diff(edges)/2)';
S.aligned_mean    = nan(opt.nBins,1);
S.misaligned_mean = nan(opt.nBins,1);
S.contrast        = nan(opt.nBins,1);
S.contrast_sem    = nan(opt.nBins,1);
S.nAligned        = zeros(opt.nBins,1);
S.nMisaligned     = zeros(opt.nBins,1);

for b = 1:opt.nBins
    inBin = (bin == b);
    a  = inBin & T.aligned == 1;
    m  = inBin & T.aligned == 0;
    S.aligned_mean(b)    = mean(T.roi_response(a), 'omitnan');
    S.misaligned_mean(b) = mean(T.roi_response(m), 'omitnan');
    S.contrast(b)        = S.aligned_mean(b) - S.misaligned_mean(b);
    S.nAligned(b)        = sum(a);
    S.nMisaligned(b)     = sum(m);

    % between-run SEM of the per-run contrast in this bin
    perRun = nan(numel(uRuns),1);
    for k = 1:numel(uRuns)
        rr = inBin & (runKey == uRuns(k));
        ra = mean(T.roi_response(rr & T.aligned==1), 'omitnan');
        rm = mean(T.roi_response(rr & T.aligned==0), 'omitnan');
        perRun(k) = ra - rm;
    end
    perRun = perRun(~isnan(perRun));
    if numel(perRun) > 1
        S.contrast_sem(b) = std(perRun) / sqrt(numel(perRun));
    end
end

% ---- report
fprintf('\nWithin-run summary (%s events, %d runs, %d bins)\n', ...
    opt.useSet, numel(uRuns), opt.nBins);
fprintf('%-10s %10s %12s %12s %8s %8s\n', ...
    'bin(frac)','aligned','misaligned','contrast','nAli','nMis');
for b = 1:opt.nBins
    fprintf('%-10.2f %10.3f %12.3f %12.3f %8d %8d\n', ...
        S.binCentre(b), S.aligned_mean(b), S.misaligned_mean(b), ...
        S.contrast(b), S.nAligned(b), S.nMisaligned(b));
end

% ---- plots
if opt.plot
    figure('Name','Within-run grid response','Color','w');

    subplot(1,2,1); hold on;
    plot(S.binCentre, S.aligned_mean,    '-o', 'LineWidth',1.5, 'DisplayName','aligned');
    plot(S.binCentre, S.misaligned_mean, '-s', 'LineWidth',1.5, 'DisplayName','misaligned');
    xlabel('within-run time (fraction)'); ylabel('ROI response (PSC)');
    title('Aligned vs misaligned over the run'); legend('Location','best'); grid on;

    subplot(1,2,2); hold on;
    errorbar(S.binCentre, S.contrast, S.contrast_sem, '-o', 'LineWidth',1.5);
    yline(0, ':');
    xlabel('within-run time (fraction)'); ylabel('aligned - misaligned (PSC)');
    title('Within-run contrast (+/- between-run SEM)'); grid on;
end
end
