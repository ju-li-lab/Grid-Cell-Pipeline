function T = check_grid_consistency(glmRoot, varargin)
% CHECK_GRID_CONSISTENCY  Does a voxel report the same grid orientation twice?
%
%   check_grid_consistency(glmRoot)
%   check_grid_consistency(glmRoot, 'Trial', 'sep4_fwhm5', 'CSV', 'grid_consistency.csv')
%
% A preprocessing-quality measure, deliberately separate from the GLM2
% contrast. If registration or unwarping is off, a voxel samples different
% tissue in run 1 and run 2, so its grid orientation stops reproducing. That
% is a preprocessing question. Whether aligned beats misaligned is an
% analysis question. This tool only answers the first.
%
% INDEPENDENT OF YOUR GLM2 SETTINGS
% ---------------------------------
% It reads ONLY the per-run images that GLM1 writes:
%
%     voxelwiseGridOri_<event>_run<N>_deg.nii
%
% and never voxelwiseGridOri_<event>_allRunsAvg_deg.nii. That matters twice
% over. First, allRunsAvg contains the run it would be compared against, so
% any agreement with it is circular - in one dataset the circular pairs read
% 78-88% where the honest run-vs-run pairs read 49-56%. Second, it means
% nothing here depends on GLM2_averageMeanGridOriAcrossRuns, so you can
% change that setting freely and these numbers stay comparable.
%
% WHAT IT REPORTS, per subject x session x hemisphere x run-pair
% --------------------------------------------------------------
%   n_mask     voxels in the ROI mask
%   n_valid    of those, voxels with an orientation in BOTH runs
%   pct_stable GridCAT's own measure: the share of voxels whose orientations
%              agree within StabilityDeg. Directly comparable with the
%              grid_metrics files. Its CHANCE LEVEL is printed - for 6-fold
%              with the default 15 degrees it is exactly 50%.
%   r_diff     resultant length of the orientation DIFFERENCE, 0-1. The
%              threshold-free version of pct_stable, and more sensitive
%              because it uses every voxel rather than counting hits.
%   z_diff     n_valid * r_diff^2. Rayleigh's statistic, so ~1 is chance and
%              >3 is roughly p<0.05. Use this to compare trials: unlike
%              pct_stable it does not throw away information, and unlike
%              r_diff it is not biased by small n.
%   r_run<N>   between-voxel coherence within each run, 0-1. This is the
%              n-independent form of the Rayleigh z in the grid_metrics
%              file, where z = n*r^2 - so r = sqrt(z/n).
%
% COMPARING PREPROCESSING TRIALS
% ------------------------------
% A mask mean changes when the voxel set changes, so a difference between
% trials can be sampling rather than signal. Two defences:
%   - n_valid is reported, so you can see when sets differ at all
%   - 'CommonMaskDir' forces every trial onto one fixed voxel set. Run each
%     trial with 'WriteValidMask', intersect the results (fslmaths -mas),
%     then re-run everything pointing at the intersection.
%
% OPTIONS
%   'Trial'         label written into the CSV so trials stack up
%   'Sessions'      default {'ses-01'}. Only session 1 by default, because an
%                   intervention between sessions could move stability for
%                   real reasons and confound the preprocessing comparison.
%                   Pass {} for all sessions once you want that.
%   'Subjects'      {} = all
%   'MaskDir'       ROI masks. Default <glmRoot>/../GLM_runauto/ROI
%   'CommonMaskDir' use these masks instead, for a fixed voxel set
%   'Event'         grid event name in the filenames (default 'translate')
%   'XFold'         default 6
%   'StabilityDeg'  default 15, which is GridCAT's own 360/XFold/4
%   'ShiftTest'     default true. Translates the mask along each axis and
%                   prints how consistency varies; a peak off centre is a
%                   coregistration shift, and its sign says which way. It
%                   prints the profile rather than a verdict - scanning 15
%                   shifts finds a peak in noise, so only a large, clearly
%                   peaked profile means anything.
%   'Shifts'        default -2:2 voxels
%   'WriteValidMask' write the valid-voxel mask next to each GLM1 dir
%   'CSV'           path to append to
%
% Hemispheres are always reported separately and never pooled - the two ErC
% masks differ in size and in how much signal survives, so a pooled number
% would hide exactly the asymmetry you want to see.
%
% Needs only base MATLAB (niftiread). SPM and GridCAT are not required.

p = inputParser;
addRequired( p, 'glmRoot',        @(x) ischar(x) || isstring(x));
addParameter(p, 'Trial',          '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Sessions',       {'ses-01'}, @iscell);
addParameter(p, 'Subjects',       {}, @iscell);
addParameter(p, 'MaskDir',        '', @(x) ischar(x) || isstring(x));
addParameter(p, 'CommonMaskDir',  '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Event',          'translate', @(x) ischar(x) || isstring(x));
addParameter(p, 'XFold',          6, @isnumeric);
addParameter(p, 'StabilityDeg',   [], @isnumeric);
addParameter(p, 'ShiftTest',      true, @islogical);
addParameter(p, 'Shifts',         -2:2, @isnumeric);
addParameter(p, 'WriteValidMask', false, @islogical);
addParameter(p, 'CSV',            '', @(x) ischar(x) || isstring(x));
parse(p, glmRoot, varargin{:});
o = p.Results;

glmRoot = char(glmRoot);
assert(isfolder(glmRoot), 'Not a folder: %s', glmRoot);

xf = o.XFold;
if isempty(o.StabilityDeg), o.StabilityDeg = 360/xf/4; end   % GridCAT's own
tolScaled  = o.StabilityDeg * xf;          % tolerance on the 360 deg circle
chancePct  = 100 * tolScaled / 180;        % two uniform angles, diff on [0,180]

maskDir = char(o.MaskDir);
if isempty(maskDir)
    maskDir = fullfile(fileparts(glmRoot), 'GLM_runauto', 'ROI');
end
if ~isempty(o.CommonMaskDir), maskDir = char(o.CommonMaskDir); end
assert(isfolder(maskDir), ['ROI masks not found: %s\n' ...
    'Pass ''MaskDir'' explicitly.'], maskDir);

glm1 = dir(fullfile(glmRoot, '*_GLM1'));
glm1 = glm1([glm1.isdir]);
assert(~isempty(glm1), 'No *_GLM1 directories in %s', glmRoot);

fprintf('\nCross-run grid orientation consistency\n');
fprintf('  GLM1 output : %s\n', glmRoot);
fprintf('  ROI masks   : %s\n', maskDir);
if ~isempty(o.Trial),    fprintf('  trial       : %s\n', char(o.Trial)); end
if ~isempty(o.Sessions), fprintf('  sessions    : %s\n', strjoin(o.Sessions, ', ')); end
fprintf('  %d-fold, stable within %g deg  ->  CHANCE = %.0f%%\n', ...
        xf, o.StabilityDeg, chancePct);
fprintf('\n%-11s %-8s %-6s %-8s %6s %7s %10s %7s %7s\n', ...
        'subject','session','hemi','runs','n_mask','n_valid','pct_stable','r_diff','z_diff');

rows = {};
for a = 1:numel(glm1)
    tok = regexp(glm1(a).name, '^(sub-[^_]+)_(ses-[^_]+)_GLM1$', 'tokens', 'once');
    if isempty(tok), continue; end
    sub = tok{1}; ses = tok{2};
    if ~isempty(o.Subjects) && ~ismember(sub, o.Subjects), continue; end
    if ~isempty(o.Sessions) && ~ismember(ses, o.Sessions), continue; end

    gdir = fullfile(glm1(a).folder, glm1(a).name);
    [ori, runNums] = load_run_orientations(gdir, char(o.Event));
    if numel(runNums) < 2
        fprintf('%-11s %-8s  (needs 2+ runs, found %d - skipped)\n', sub, ses, numel(runNums));
        continue;
    end

    masks = dir(fullfile(maskDir, sprintf('r%s_%s_hemi-*_mask.nii', sub, ses)));
    if isempty(masks)
        fprintf('%-11s %-8s  (no ROI masks for this session - skipped)\n', sub, ses);
        continue;
    end

    for m = 1:numel(masks)
        M = niftiread(fullfile(masks(m).folder, masks(m).name));
        if numel(M) ~= numel(ori{1})
            error(['Mask and orientation images are on different grids:\n' ...
                   '  %s  (%d voxels)\n  orientation image (%d voxels)'], ...
                   masks(m).name, numel(M), numel(ori{1}));
        end
        hemi    = hemi_of(masks(m).name);
        maskIdx = find(M(:) > 0.5);
        if isempty(maskIdx), continue; end

        % between-voxel coherence per run, on the voxels that have a value
        rRun = nan(1, numel(runNums));
        for k = 1:numel(runNums)
            v = ori{k}(maskIdx);
            v = v(~isnan(v));
            if ~isempty(v), rRun(k) = resultant(v * xf); end
        end

        for i = 1:numel(runNums)
            for j = i+1:numel(runNums)
                [st, rd, zd, nValid] = compare_runs( ...
                    ori{i}(maskIdx), ori{j}(maskIdx), xf, tolScaled);
                cmp = sprintf('%d_vs_%d', runNums(i), runNums(j));
                fprintf('%-11s %-8s %-6s %-8s %6d %7d %9.1f%% %7.3f %7.2f\n', ...
                        sub, ses, hemi, cmp, numel(maskIdx), nValid, st, rd, zd);
                rows{end+1,1} = struct('trial', char(o.Trial), 'subject', sub, ...
                    'session', ses, 'hemi', hemi, 'comparison', cmp, ...
                    'n_mask', numel(maskIdx), 'n_valid', nValid, ...
                    'pct_stable', st, 'r_diff', rd, 'z_diff', zd, ...
                    'chance_pct', chancePct, ...
                    'r_run_a', rRun(i), 'r_run_b', rRun(j)); %#ok<AGROW>
            end
        end

        if o.ShiftTest
            shift_test(ori, runNums, M, xf, tolScaled, o.Shifts, sub, ses, hemi);
        end
        if o.WriteValidMask
            write_valid_mask(gdir, masks(m), M, ori, sub, ses, hemi);
        end
    end
end

T = summarise(rows, chancePct);
if ~isempty(o.CSV), write_csv(char(o.CSV), rows); end
end


% ---------------------------------------------------------------------------

function [ori, runNums] = load_run_orientations(gdir, event)
% Per-run images ONLY. allRunsAvg is excluded on purpose: comparing a run
% against an average that contains it is circular, and depending on it would
% tie these numbers to GLM2_averageMeanGridOriAcrossRuns.
f = dir(fullfile(gdir, sprintf('voxelwiseGridOri_%s_run*_deg.nii', event)));
ori = {}; runNums = [];
for k = 1:numel(f)
    n = regexp(f(k).name, '_run(\d+)_deg\.nii$', 'tokens', 'once');
    if isempty(n), continue; end
    runNums(end+1) = str2double(n{1}); %#ok<AGROW>
    ori{end+1}     = double(niftiread(fullfile(f(k).folder, f(k).name))); %#ok<AGROW>
end
[runNums, order] = sort(runNums);
ori = ori(order);
end


function [pctStable, rDiff, zDiff, nValid] = compare_runs(a, b, xf, tolScaled)
% Orientations are stored on [0, 360/xf). Scaling by xf puts them on the full
% circle, which is where the circular statistics belong.
ok = ~isnan(a) & ~isnan(b);
nValid = sum(ok);
if nValid == 0
    [pctStable, rDiff, zDiff] = deal(NaN); return;
end
d = abs_ang_diff(a(ok) * xf, b(ok) * xf, 360);     % folded to [0,180]
pctStable = 100 * sum(d <= tolScaled) / nValid;
rDiff     = resultant((a(ok) - b(ok)) * xf);        % signed difference
zDiff     = nValid * rDiff^2;                       % Rayleigh: ~1 is chance
end


function r = resultant(deg)
% Mean resultant length of angles given in degrees. 0 = uniform, 1 = identical.
r = abs(mean(exp(1i * deg(:) * pi/180)));
end


function d = abs_ang_diff(a, b, maxAng)
n = mod(a - b, maxAng);
d = min(maxAng - n, n);
end


function shift_test(ori, runNums, M, xf, tolScaled, shifts, sub, ses, hemi)
% Translate the mask along each axis and print how consistency varies. If
% the mask and the functional data disagree about where the ROI is - a
% coregistration shift - the peak sits off centre, and its sign says which
% way.
%
% It prints the profile rather than passing a verdict, on purpose. Scanning
% five shifts on three axes WILL find a peak in pure noise (the maximum of
% 15 draws from Rayleigh's null sits around z = 3), so only a large, clearly
% peaked profile means anything. A flat profile means the scan found nothing,
% which is the answer you want.
axes_ = 'xyz';
for ax = 1:3
    zz = nan(size(shifts));
    for s = 1:numel(shifts)
        v = [0 0 0]; v(ax) = shifts(s);
        zz(s) = median_z(ori, runNums, find(circshift(M, v) > 0.5), xf, tolScaled);
    end
    % prefer no shift when the profile is flat, so ties do not read as an offset
    z0 = zz(shifts == 0);
    [mx, k] = max(zz);
    if ~isempty(z0) && z0 >= 0.98 * mx
        peak = 0;
    else
        peak = shifts(k);
    end
    fprintf('%-11s %-8s %-6s   shift %c:', sub, ses, hemi, axes_(ax));
    fprintf(' %7.2f', zz);
    fprintf('   peak %+d\n', peak);
end
end


function z = median_z(ori, runNums, idx, xf, tolScaled)
zs = [];
for i = 1:numel(runNums)
    for j = i+1:numel(runNums)
        [~, ~, zd] = compare_runs(ori{i}(idx), ori{j}(idx), xf, tolScaled);
        zs(end+1) = zd; %#ok<AGROW>
    end
end
z = median(zs, 'omitnan');
end


function write_valid_mask(gdir, maskFile, M, ori, sub, ses, hemi)
% The voxels that actually have an orientation in every run. Intersect these
% across trials to force one fixed voxel set.
V = M > 0.5;
for k = 1:numel(ori)
    V = V & ~isnan(ori{k});
end
info = niftiinfo(fullfile(maskFile.folder, maskFile.name));
info.Datatype = 'uint8';
out = fullfile(gdir, sprintf('validvox_%s_%s_%s.nii', sub, ses, hemi));
niftiwrite(uint8(V), out, info);
end


function h = hemi_of(name)
n = lower(name);
if     contains(n, 'bilat'), h = 'bilat';
elseif contains(n, 'left'),  h = 'left';
elseif contains(n, 'right'), h = 'right';
else,                        h = 'other';
end
end


function T = summarise(rows, chancePct)
if isempty(rows)
    fprintf('\nNothing measured.\n'); T = table(); return;
end
hemi = cellfun(@(r) string(r.hemi),   rows);
st   = cellfun(@(r) r.pct_stable,     rows);
zd   = cellfun(@(r) r.z_diff,         rows);

fprintf('\n%d run-pairs measured.  Chance for pct_stable is %.0f%%.\n\n', numel(rows), chancePct);
fprintf('%-8s %8s %12s %12s %s\n', 'hemi', 'n pairs', 'median %', 'median z', 'verdict');
for h = unique(hemi)'
    k = hemi == h;
    mz = median(zd(k), 'omitnan');
    if     mz > 3,   v = 'above chance';
    elseif mz > 1.5, v = 'marginal';
    else,            v = 'AT CHANCE - orientations do not reproduce';
    end
    fprintf('%-8s %8d %11.1f%% %12.2f %s\n', h, sum(k), median(st(k),'omitnan'), mz, v);
end
fprintf(['\nz_diff is Rayleigh''s statistic on the orientation difference:\n' ...
         '~1 is chance, >3 is roughly p<0.05. Compare the SAME subject across\n' ...
         'preprocessing trials - a subject with no grid signal sits at chance\n' ...
         'however good the preprocessing, so this is not a between-subject measure.\n']);

h  = unique(hemi);
np = arrayfun(@(x) sum(hemi == x), h);
mp = arrayfun(@(x) median(st(hemi == x), 'omitnan'), h);
mz = arrayfun(@(x) median(zd(hemi == x), 'omitnan'), h);
T  = table(h, np, mp, mz, 'VariableNames', ...
           {'hemi','n_pairs','median_pct_stable','median_z_diff'});
end


function write_csv(path, rows)
newFile = ~isfile(path);
fid = fopen(path, 'a');
assert(fid > 0, 'Cannot write %s', path);
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
if newFile
    fprintf(fid, ['trial,subject,session,hemi,comparison,n_mask,n_valid,' ...
                  'pct_stable,chance_pct,r_diff,z_diff,r_run_a,r_run_b\n']);
end
for i = 1:numel(rows)
    r = rows{i};
    fprintf(fid, '%s,%s,%s,%s,%s,%d,%d,%.2f,%.2f,%.4f,%.4f,%.4f,%.4f\n', ...
            r.trial, r.subject, r.session, r.hemi, r.comparison, ...
            r.n_mask, r.n_valid, r.pct_stable, r.chance_pct, ...
            r.r_diff, r.z_diff, r.r_run_a, r.r_run_b);
end
fprintf('\nAppended %d rows to %s\n', numel(rows), path);
end
