function roi_tsnr_filter(derivRoot, varargin)
%ROI_TSNR_FILTER  Drop low-tSNR outlier voxels from each ROI mask.
%
%   roi_tsnr_filter(derivRoot, 'RoiDir', D)
%   roi_tsnr_filter(derivRoot, 'RoiDir', D, 'DryRun', true)
%
% Exclusion rule: within each subject x session, pool the tSNR values of all
% ROI voxels (left + right) and drop every voxel more than 1.5 interquartile
% ranges below the first quartile (tSNR < Q1 - 1.5*IQR). A standard outlier
% fence, computed within each subject, so no threshold is carried between
% subjects. Runs on every session.
%
% Writes filtered copies of the left/right masks (filenames unchanged),
% rebuilds bilat as the union of the two filtered sides, and writes the tSNR
% image it used. Default output: <RoiDir>_tsnr, which is what 03_gridcat.sh
% reads when ROI_TSNR_FILTER=1 in config.sh.
%
% tSNR here is mean(raw) / SD(residual after a DCT high-pass at the GLM's own
% cutoff) - the noise the GLM actually sees. It is direction-blind, so
% selecting on it cannot bias the GLM2 contrast. Use UNSMOOTHED functionals
% (prefix 'u'/'ua'): smoothed tSNR measures distance to tissue edge, not the
% voxel.
%
% Options:
%   'RoiDir'    flat dir of r<sub>_<ses>_hemi-*_label-*_mask.nii.
%               Default <derivRoot>/../GLM_runauto/ROI
%   'OutDir'    default <RoiDir>_tsnr
%   'Subjects'  {} = all
%   'Sessions'  {} = all (the default)
%   'Prefix'    functional prefix, default 'u'. Never a smoothed series.
%   'Tasks'     default [1 2 3], matching TASKS in config.sh
%   'Combine'   'mean' (default) or 'min' across runs
%   'HPFSec'    high-pass seconds, default 128, matching GLM.HPF_cutoff
%   'TR'        default 2.06; read from the header when it looks sane
%   'MinVoxels' refuse to write a mask below this, default 5
%   'DryRun'    report only, write nothing
%   'CSV'       default <OutDir>/roi_tsnr.csv
%
% One subject x session per call is the unit of parallelism - see
% run_roi_tsnr_filter.sh. Nothing in the pipeline calls this automatically.

p = inputParser;
p.addParameter('RoiDir',    '');
p.addParameter('OutDir',    '');
p.addParameter('Subjects',  {});
p.addParameter('Sessions',  {});
p.addParameter('Prefix',    'u');
p.addParameter('Tasks',     [1 2 3]);
p.addParameter('Combine',   'mean');
p.addParameter('HPFSec',    128);
p.addParameter('TR',        2.06);
p.addParameter('MinVoxels', 5);
p.addParameter('DryRun',    false);
p.addParameter('Trial',     '');
p.addParameter('CSV',       '');
p.parse(varargin{:});
o = p.Results;

if ~isfolder(derivRoot), error('No such derivatives root: %s', derivRoot); end
if ~ismember(lower(o.Combine), {'mean', 'min'})
    error('Combine must be ''mean'' or ''min'' (got: %s)', o.Combine);
end

roiDir = o.RoiDir;
if isempty(roiDir)
    roiDir = fullfile(fileparts(derivRoot), 'GLM_runauto', 'ROI');
end
if ~isfolder(roiDir), error('No such ROI directory: %s', roiDir); end

outDir = o.OutDir;
if isempty(outDir), outDir = [roiDir '_tsnr']; end
if ~o.DryRun && ~isfolder(outDir), mkdir(outDir); end
csvPath = o.CSV;
if isempty(csvPath), csvPath = fullfile(outDir, 'roi_tsnr.csv'); end

% ---- discover subject x session from the ROI masks --------------------
m = dir(fullfile(roiDir, 'r*_hemi-left_label-*_mask.nii'));
keys = {};
for a = 1:numel(m)
    tok = regexp(m(a).name, '^r(sub-[^_]+)_(ses-[^_]+)_hemi-left_label-([^_]+)_mask\.nii$', ...
                 'tokens', 'once');
    if isempty(tok), continue; end
    if ~isempty(o.Subjects) && ~ismember(tok{1}, o.Subjects), continue; end
    if ~isempty(o.Sessions) && ~ismember(tok{2}, o.Sessions), continue; end
    keys = [keys; tok(:)'];                                        %#ok<AGROW>
end
if isempty(keys)
    error(['No r<sub>_<ses>_hemi-left_label-*_mask.nii under %s ' ...
           'matching the requested subjects/sessions.'], roiDir);
end

fprintf('\n============================================\n');
fprintf('  ROI tSNR filter\n');
fprintf('  Functionals : %s  (%ssub-*, unsmoothed)\n', derivRoot, o.Prefix);
fprintf('  ROI masks   : %s\n', roiDir);
fprintf('  Rule        : drop voxels with tSNR < Q1 - 1.5*IQR, per subject x session\n');
fprintf('  Runs        : tasks [%s], combined by %s\n', num2str(o.Tasks), lower(o.Combine));
fprintf('  High-pass   : DCT, %g s\n', o.HPFSec);
if o.DryRun
    fprintf('  DRY RUN - reporting only, nothing will be written.\n');
else
    fprintf('  Out         : %s\n', outDir);
end
fprintf('============================================\n\n');

rows = {};
for k = 1:size(keys, 1)
    sub = keys{k,1}; ses = keys{k,2}; label = keys{k,3};

    funcs = find_runs(derivRoot, sub, ses, o.Prefix, o.Tasks);
    if isempty(funcs)
        fprintf('%-11s %-8s  no %ssub-* runs found - skipped\n', sub, ses, o.Prefix);
        continue;
    end

    refPath = fullfile(roiDir, sprintf('r%s_%s_hemi-left_label-%s_mask.nii', sub, ses, label));
    info    = niftiinfo(refPath);

    [tsnr, funcInfo, brainTsnr] = tsnr_volume(funcs, o.HPFSec, o.TR, lower(o.Combine));
    check_same_grid(info, funcInfo, refPath, funcs{1});

    % ---- pooled distribution over both hemispheres -> outlier fence ----
    sides  = {'left', 'right'};
    roiIdx = struct('left', [], 'right', []);
    pooled = [];
    for h = 1:2
        mPath = fullfile(roiDir, sprintf('r%s_%s_hemi-%s_label-%s_mask.nii', ...
                                         sub, ses, sides{h}, label));
        if ~isfile(mPath), continue; end
        M = niftiread(mPath);
        if numel(M) ~= numel(tsnr)
            error(['ROI mask and the functionals are on different grids:\n' ...
                   '  %s (%d voxels)\n  functional (%d voxels)'], ...
                   mPath, numel(M), numel(tsnr));
        end
        roiIdx.(sides{h}) = find(M(:) > 0.5);
        tv = tsnr(roiIdx.(sides{h}));
        pooled = [pooled; tv(isfinite(tv))];                       %#ok<AGROW>
    end
    if isempty(pooled)
        fprintf('%-11s %-8s  no finite tSNR values in the ROI - skipped\n', sub, ses);
        continue;
    end

    q1  = quartile(pooled, 0.25);
    q3  = quartile(pooled, 0.75);
    cut = q1 - 1.5 * (q3 - q1);

    fprintf('%s %s  label-%s  %d runs  brain tSNR %.1f  Q1 %.1f  IQR %.1f  cut %.1f\n', ...
            sub, ses, label, numel(funcs), brainTsnr, q1, q3 - q1, cut);

    kept = struct('left', [], 'right', []);
    for h = 1:2
        side = sides{h};
        idx  = roiIdx.(side);
        if isempty(idx)
            fprintf('   %-5s  no mask - skipped\n', side);
            continue;
        end
        tv  = tsnr(idx);
        sel = isfinite(tv) & tv >= cut;

        fprintf('   %-5s  %3d voxels, tSNR median %.1f -> dropping %d outlier(s)\n', ...
                side, numel(idx), median(tv(isfinite(tv))), sum(~sel));

        rows = [rows; {o.Trial, sub, ses, side, label, numel(funcs), ...
                       numel(idx), sum(sel), 100 * sum(sel) / numel(idx), cut, ...
                       median(tv(sel), 'omitnan'), ...
                       median(tv(isfinite(tv) & ~sel), 'omitnan'), ...
                       median(tv(isfinite(tv)), 'omitnan'), brainTsnr}];  %#ok<AGROW>

        if sum(sel) < o.MinVoxels
            warning('roi_tsnr_filter:tooFew', ...
                '%s %s hemi-%s would keep only %d voxels (MinVoxels %d). NOT WRITTEN.', ...
                 sub, ses, side, sum(sel), o.MinVoxels);
            continue;
        end
        kept.(side) = idx(sel);

        if ~o.DryRun
            write_mask(info, size(tsnr), idx(sel), ...
                fullfile(outDir, sprintf('r%s_%s_hemi-%s_label-%s_mask', ...
                                         sub, ses, side, label)));
        end
    end

    if ~o.DryRun
        % bilat is rebuilt from the FILTERED sides, never copied
        both = [kept.left; kept.right];
        if numel(both) >= o.MinVoxels
            write_mask(info, size(tsnr), both, ...
                fullfile(outDir, sprintf('r%s_%s_hemi-bilat_label-%s_mask', ...
                                         sub, ses, label)));
        end
        write_vol(info, tsnr, fullfile(outDir, sprintf('%s_%s_tsnr', sub, ses)));
    end
end

if ~o.DryRun
    append_csv(csvPath, rows, ['trial,subject,session,hemi,label,n_runs,' ...
        'n_roi,n_kept,pct_kept,cut_tsnr,med_kept,med_dropped,med_all,brain_tsnr'], ...
        '%s,%s,%s,%s,%s,%d,%d,%d,%.2f,%.3f,%.3f,%.3f,%.3f,%.3f\n');
    write_info(outDir, derivRoot, roiDir, o);
    fprintf('\nmasks + tSNR images -> %s\n', outDir);
    fprintf('%d rows -> %s\n', size(rows, 1), csvPath);
    fprintf('\nTo use them: set ROI_TSNR_FILTER=1 in config.sh.\n');
else
    fprintf('\nDry run - nothing written.\n');
end

end

% ======================================================================

function q = quartile(x, p)
% Linear-interpolation quantile (R type 7), no Statistics Toolbox needed.
x = sort(x(:));
n = numel(x);
if n == 1, q = x; return; end
h = (n - 1) * p + 1;
lo = floor(h);
hi = ceil(h);
q = x(lo) + (h - lo) * (x(hi) - x(lo));
end

function funcs = find_runs(derivRoot, sub, ses, prefix, tasks)
% One 4D file per task. Anchored on the prefix so 'u' cannot pick up 'ua'/'su'.
funcs = {};
fdir = fullfile(derivRoot, sub, ses, 'func');
for t = tasks(:)'
    g = dir(fullfile(fdir, sprintf('%ssub-*_%s_task-%d*_bold.nii', prefix, ses, t)));
    g = g(~startsWith({g.name}, 'mean'));
    for a = 1:numel(g)
        nm = g(a).name;
        if ~strcmp(nm(1:numel(prefix)), prefix) || ...
           ~strncmp(nm(numel(prefix)+1:end), 'sub-', 4)
            continue;
        end
        funcs{end+1} = fullfile(g(a).folder, nm);                  %#ok<AGROW>
    end
end
end

function [tsnr, funcInfo, brainTsnr] = tsnr_volume(funcs, hpfSec, TR, how)
% mean(raw) / SD(residual after a DCT high-pass), per run, then combined.
CHUNK = 20000;
acc = []; info = []; nUsed = 0;
for r = 1:numel(funcs)
    if isempty(info), info = niftiinfo(funcs{r}); end
    V  = niftiread(funcs{r});
    sz = size(V);
    if numel(sz) < 4 || sz(4) < 2
        warning('roi_tsnr_filter:not4D', 'Not a 4D run, skipped: %s', funcs{r});
        continue;
    end
    n  = sz(4);
    tr = TR;
    pd = info.PixelDimensions;
    if numel(pd) >= 4 && pd(4) > 0.1 && pd(4) < 20, tr = double(pd(4)); end

    Y  = reshape(single(V), [], n);                 clear V
    mu = mean(Y, 2);
    brain = mu > mean(mu) / 8;

    X0 = dct_basis(n, fix(2 * (n * tr) / hpfSec + 1));
    df = max(n - size(X0, 2), 1);

    sd  = zeros(size(mu), 'single');
    idx = find(brain);
    for c = 1:CHUNK:numel(idx)
        j  = idx(c:min(c + CHUNK - 1, numel(idx)));
        Yc = double(Y(j, :))';                      % nScan x nVox
        Rc = Yc - X0 * (X0' * Yc);                  % X0 is orthonormal
        sd(j) = single(sqrt(sum(Rc.^2, 1) / df))';
    end
    clear Y

    t = nan(size(mu), 'single');
    good = brain & sd > 0;
    t(good) = mu(good) ./ sd(good);
    t = reshape(t, sz(1:3));

    nUsed = nUsed + 1;
    if isempty(acc)
        acc = t;
    elseif strcmp(how, 'min')
        acc = min(acc, t);
    else
        acc = acc + t;
    end
end
if isempty(acc), error('No usable 4D run among the %d matched files.', numel(funcs)); end
if strcmp(how, 'mean'), acc = acc / nUsed; end
tsnr = acc;
funcInfo = info;
brainTsnr = median(tsnr(isfinite(tsnr)));
end

function check_same_grid(maskInfo, funcInfo, maskPath, funcPath)
a = maskInfo.Transform.T; b = funcInfo.Transform.T;
if ~isequal(maskInfo.ImageSize(1:3), funcInfo.ImageSize(1:3)) || ~all(abs(a(:) - b(:)) < 1e-3)
    error(['Mask and functional are not on the same grid:\n  %s\n  %s\n' ...
           'Use the resliced r*_mask.nii from the same derivatives tree as ' ...
           'the functionals.'], maskPath, funcPath);
end
end

function C = dct_basis(n, k)
% Orthonormal DCT-II basis, column 1 constant - the set SPM's high-pass removes.
k = max(min(k, n), 1);
C = zeros(n, k);
C(:, 1) = 1 / sqrt(n);
j = 2:k;
C(:, j) = sqrt(2 / n) * cos(pi * (2 * (0:n-1)' + 1) * (j - 1) / (2 * n));
end

function write_mask(refInfo, sz, idx, stemNoExt)
M = zeros(sz, 'single');
M(idx) = 1;
write_vol(refInfo, M, stemNoExt);
end

function write_vol(refInfo, X, stemNoExt)
info = refInfo;
info.ImageSize             = size(X);
info.PixelDimensions       = refInfo.PixelDimensions(1:3);
info.Datatype              = 'single';
info.raw.datatype          = 16;
info.BitsPerPixel          = 32;
info.Description           = 'roi_tsnr_filter';
info.MultiplicativeScaling = 1;
info.AdditiveOffset        = 0;
niftiwrite(single(X), stemNoExt, info);
end

function append_csv(path, rows, hdr, fmt)
if isempty(rows), return; end
newFile = ~isfile(path);
fid = fopen(path, 'a');
if fid < 0, error('Cannot write CSV: %s', path); end
if newFile, fprintf(fid, '%s\n', hdr); end
for i = 1:size(rows, 1), fprintf(fid, fmt, rows{i,:}); end
fclose(fid);
end

function write_info(outDir, derivRoot, roiDir, o)
fid = fopen(fullfile(outDir, 'filter_info.txt'), 'w');
if fid < 0, return; end
fprintf(fid, 'roi_tsnr_filter\n');
fprintf(fid, 'written          %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, 'functionals      %s (prefix %ssub-*, unsmoothed)\n', derivRoot, o.Prefix);
fprintf(fid, 'source ROI dir   %s\n', roiDir);
fprintf(fid, 'rule             drop voxels with tSNR < Q1 - 1.5*IQR (pooled L+R, per subject x session)\n');
fprintf(fid, 'tasks            [%s]\n', num2str(o.Tasks));
fprintf(fid, 'runs combined by %s\n', lower(o.Combine));
fprintf(fid, 'high-pass        DCT, %g s\n', o.HPFSec);
fprintf(fid, 'min voxels       %d\n', o.MinVoxels);
fprintf(fid, 'trial            %s\n', o.Trial);
fprintf(fid, '\nbilat masks are the union of the FILTERED left and right.\n');
fclose(fid);
end
