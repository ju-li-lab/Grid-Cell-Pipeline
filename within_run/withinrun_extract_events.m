function T = withinrun_extract_events(glm1Dir, roiMask, outCsv, varargin)
% WITHINRUN_EXTRACT_EVENTS  Per-event ROI response for looking *inside* a run.
%
% GridCAT collapses each run into a single aligned-vs-misaligned contrast
% value. That value is one GLM beta contrast estimated from ~19 aligned and
% ~19 misaligned events pooled together, so by construction it cannot tell
% you how the grid signal behaves early vs. late in a run. This function
% produces the raw substrate you need to see that: ONE response number per
% grid event, tagged with its within-run time and its alignment to the mean
% grid orientation. Everything downstream (early/late split, sliding window,
% regression on time) is then just filtering/plotting this table -- no new
% GLM required.
%
% WHAT IT DOES (all transparent, no SPM GLM re-fit):
%   1. Reproduces the SAME mean grid orientation phi the pipeline uses,
%      by calling GridCAT's own calcMeanGridOri on the stored GLM1 result.
%   2. Reads each run's event table (GridCAT's own readEventTable).
%   3. Extracts the ROI-averaged BOLD timeseries directly from the
%      preprocessed functional volumes.
%   4. For each grid event, samples the ROI response in a post-onset window
%      (a crude peristimulus read-out -- see LIMITATIONS below).
%   5. Writes a tidy long-format CSV: one row per event.
%
% This is deliberately the simplest fundamental mechanism ("look at the
% signal around each event"), not the most rigorous one. See README.md in
% this folder for the ladder toward single-trial (LSS) and temporal-bin
% GLMs once you know what you are looking for.
%
% USAGE
%   % After a normal pipeline run, each session has a *_GLM1 folder whose
%   % GridCAT_GLM1.mat stores cfg (functional scans, event tables, TR...)
%   % and the estimated voxelwise orientations.
%   T = withinrun_extract_events( ...
%         '/.../GLM_output/sub-01s13_ses-01_GLM1', ...   % GLM1 result dir
%         '/.../ROI/rsub-01s13_ses-01_ErC-bilat.nii', ...% same mask as GLM2
%         '/.../within_run/sub-01s13_ses-01_bilat_events.csv');
%
% OPTIONS (name/value) -- defaults chosen to match pipeline_config.cfg
%   'useWeighting'      1     amplitude-weight voxels in phi (USE_VOXEL_WEIGHTING)
%   'avgOriAcrossRuns'  true  phi from allRunsAvg (AVG_ORI_ACROSS_RUNS=1);
%                             false -> per-run phi
%   'eventUsageGLM1'    4     specifier used to tag GLM1 (estimation) events
%   'eventUsageGLM2'    5     specifier used to tag GLM2 (test) events
%   'gridEventName'     ''    '' = all grid event types; else e.g. 'translate'
%   'respWindow'        [4 8] seconds post-onset to average as the response
%   'baseWindow'        []    e.g. [-2 0] to subtract a pre-onset baseline;
%                             [] = no baseline subtraction
%   'detrendOrder'      2     polynomial order removed from ROI timeseries
%                             (light-touch stand-in for the GLM high-pass)
%
% OUTPUT columns
%   sub, ses, run, event_in_run, onset_s, within_run_frac,
%   orientation_deg, mean_grid_ori_deg, alignment, aligned, in_glm1, in_glm2,
%   roi_response
%     alignment = cos(xFold*(theta - phi)); aligned = alignment >= 0
%     (identical definition to GridCAT's aligned/misaligned split).
%
% LIMITATIONS (important -- this is a first-look tool):
%   * The response is read straight off the ROI-mean timeseries. With ~8 s
%     between events the HRFs overlap, so a single event's number is
%     contaminated by its neighbours. It is fine for spotting run-level
%     trends (the neighbour contamination averages out across many events),
%     but it is NOT a clean single-trial estimate.
%   * Motion and non-grid task events are not regressed out here. The GridCAT
%     contrast removes them; this does not. Treat absolute values with care;
%     the aligned-minus-misaligned *difference* cancels most shared nuisance.
%   * For a rigorous per-event estimate, move to an LSS beta-series GLM
%     (see README.md). This tool exists to let you look first.

% ---- reuse the GridCAT toolbox functions (readEventTable, calcMeanGridOri,
%      loadImage_SPM) and SPM (spm_vol, spm_get_data). They must be on the
%      path, exactly as for a normal pipeline run.
assert(exist('calcMeanGridOri','file')==2, 'GridCAT not on path (calcMeanGridOri missing).');
assert(exist('spm_vol','file')==2,        'SPM not on path (spm_vol missing).');

% ---- parse options
p = inputParser;
p.addParameter('useWeighting', 1);
p.addParameter('avgOriAcrossRuns', true);
p.addParameter('eventUsageGLM1', 4);
p.addParameter('eventUsageGLM2', 5);
p.addParameter('gridEventName', '');
p.addParameter('respWindow', [4 8]);
p.addParameter('baseWindow', []);
p.addParameter('detrendOrder', 2);
p.parse(varargin{:});
opt = p.Results;

% ---- load the stored GLM1 result (cfg + voxelwise orientations)
glm1File = fullfile(glm1Dir, 'GridCAT_GLM1.mat');
assert(exist(glm1File,'file')==2, 'GridCAT_GLM1.mat not found in %s', glm1Dir);
G = load(glm1File);                       % fields: cfg, gridEventType
cfg = G.cfg;
xFold = cfg.GLM.xFoldSymmetry;
TR    = cfg.GLM.TR;

% ---- derive sub/ses labels from the first event table name (for the CSV)
[sub, ses] = parse_sub_ses(cfg.rawData.run(1).eventTable_file);

% ---- load ROI mask once and turn it into voxel coordinates for spm_get_data
[maskHdr, maskVol] = loadImage_SPM(roiMask);
roiLin = find(maskVol > 0);
[ix, iy, iz] = ind2sub(size(maskVol), roiLin);
roiXYZ = [ix, iy, iz]';                    % 3 x nVox voxel coordinates
fprintf('ROI %s: %d voxels\n', roiMask, numel(roiLin));

rows = {};   % accumulate output rows

nRun = numel(cfg.rawData.run);
for r = 1:nRun

    % -------- read this run's events (GridCAT's own reader) --------
    ev = readEventTable(cfg.rawData.run(r).eventTable_file);

    % which grid event types to process
    if isempty(opt.gridEventName)
        typeIdx = 1:numel(ev.gridEventType);
    else
        typeIdx = find(strcmpi({ev.gridEventType.eventName}, opt.gridEventName));
        assert(~isempty(typeIdx), 'gridEventName "%s" not found in run %d', opt.gridEventName, r);
    end

    % -------- load the ROI-mean BOLD timeseries for this run --------
    V = spm_vol(char(cfg.rawData.run(r).functionalScans));   % nVol spm_vol structs
    assert(isequal(V(1).dim(:)', maskHdr.dim(:)'), ...
        'Run %d functional dim %s != mask dim %s. Reslice the mask to functional space.', ...
        r, mat2str(V(1).dim(:)'), mat2str(maskHdr.dim(:)'));
    Y = spm_get_data(V, roiXYZ);           % nVol x nVox
    roiTS = mean(Y, 2, 'omitnan');         % nVol x 1  ROI-mean signal
    nVol  = numel(roiTS);
    runDur = nVol * TR;

    % light detrend + convert to percent signal change (see LIMITATIONS)
    roiPSC = detrend_psc(roiTS, opt.detrendOrder);

    % volume centre times (s from run start)
    volT = ((0:nVol-1) + 0.5) * TR;

    % -------- per grid event type --------
    for t = typeIdx
        eventName = ev.gridEventType(t).eventName;
        thetaDeg  = ev.gridEventType(t).orientations(:);     % nEvents x 1, degrees
        nEvents   = numel(thetaDeg);

        % onsets for this event type (readEventTable stores onsets per unique
        % name in ev.onsets; match by name)
        onsetIdx  = find(strcmp(ev.names, eventName), 1);
        onsets    = ev.onsets{onsetIdx}(:);                  % nEvents x 1, seconds
        assert(numel(onsets) == nEvents, 'onset/orientation count mismatch, run %d', r);

        % same phi the pipeline used: find this event type in the stored GLM1
        gi = find(strcmp(eventName, {G.gridEventType.eventName}), 1);
        assert(~isempty(gi), 'event "%s" not in GLM1 result', eventName);
        if opt.avgOriAcrossRuns
            oriSrc = G.gridEventType(gi).allRunsAvg;
        else
            oriSrc = G.gridEventType(gi).run(r);
        end
        phi = calcMeanGridOri(xFold, roiMask, oriSrc, opt.useWeighting);   % radians

        % alignment, identical definition to generateMultiCondFile
        alignment = cos(xFold * (deg2rad(thetaDeg) - phi));   % nEvents x 1
        aligned   = double(alignment >= 0);

        % GLM1 / GLM2 membership tags (replicate specifyEventUsageForGLM)
        inGLM1 = usage_logical(opt.eventUsageGLM1, nEvents, r);
        inGLM2 = usage_logical(opt.eventUsageGLM2, nEvents, r);

        % -------- per-event ROI response (peristimulus window mean) --------
        for e = 1:nEvents
            resp = window_mean(roiPSC, volT, onsets(e) + opt.respWindow);
            if ~isempty(opt.baseWindow)
                resp = resp - window_mean(roiPSC, volT, onsets(e) + opt.baseWindow);
            end
            rows(end+1, :) = { sub, ses, r, e, onsets(e), onsets(e)/runDur, ...
                thetaDeg(e), rad2deg(phi), alignment(e), aligned(e), ...
                inGLM1(e), inGLM2(e), resp }; %#ok<AGROW>
        end
    end
end

% ---- assemble table and write CSV
T = cell2table(rows, 'VariableNames', { ...
    'sub','ses','run','event_in_run','onset_s','within_run_frac', ...
    'orientation_deg','mean_grid_ori_deg','alignment','aligned', ...
    'in_glm1','in_glm2','roi_response'});

if nargin >= 3 && ~isempty(outCsv)
    outDir = fileparts(outCsv);
    if ~isempty(outDir) && ~exist(outDir,'dir'); mkdir(outDir); end
    writetable(T, outCsv);
    fprintf('Wrote %d event rows -> %s\n', height(T), outCsv);
end
end


% ===================== local helpers =====================

function [sub, ses] = parse_sub_ses(eventFile)
[~, name] = fileparts(eventFile);
tok = regexp(name, '(sub-[^_]+)_(ses-[^_]+)', 'tokens', 'once');
if isempty(tok); sub = 'sub-?'; ses = 'ses-?'; else; sub = tok{1}; ses = tok{2}; end
end

function y = detrend_psc(ts, order)
% Remove a low-order polynomial drift (stand-in for the GLM high-pass) and
% express fluctuations as percent signal change relative to the run mean.
ts = ts(:);
n  = numel(ts);
m  = mean(ts, 'omitnan');
if order >= 1
    x = linspace(-1, 1, n)';
    X = zeros(n, order);
    for pw = 1:order; X(:,pw) = x.^pw; end   % no constant column -> mean preserved
    good = ~isnan(ts);
    beta = X(good,:) \ (ts(good) - m);
    ts = ts - X*beta;
end
y = 100 * (ts - m) / m;                        % percent signal change, centred at 0
end

function v = window_mean(psc, volT, win)
% Mean of the ROI signal over volumes whose centre time falls in [win(1) win(2)].
sel = volT >= win(1) & volT < win(2);
if any(sel); v = mean(psc(sel), 'omitnan'); else; v = NaN; end
end

function ids = usage_logical(spec, nEvents, runIdx)
% Replicates GridCAT/specifyEventUsageForGLM so the tags here match exactly
% which events the pipeline fed to GLM1 / GLM2.
ids = zeros(nEvents, 1);
switch spec
    case 2, ids(1:round(nEvents/2)) = 1;                        % first half
    case 3, ids(:) = 1; ids(1:round(nEvents/2)) = 0;           % second half
    case 4, ids(1:2:nEvents) = 1;                              % odd events
    case 5, ids(2:2:nEvents) = 1;                              % even events
    case 6, if mod(runIdx,2); ids(:) = 1; end                 % odd runs
    case 7, if ~mod(runIdx,2); ids(:) = 1; end                % even runs
    case 8, ids(:) = 1;                                        % all events
    otherwise, ids(:) = NaN;   % 9 (from-table) not replicated here
end
end
