function roi_tsnr_filter(derivRoot, varargin)
%ROI_TSNR_FILTER  Keep only the high-tSNR voxels of each ROI mask.
%
%   roi_tsnr_filter(derivRoot, 'RoiDir', D)                    % report only
%   roi_tsnr_filter(derivRoot, 'RoiDir', D, 'KeepPct', [60 40])  % L 60%, R 40%
%   roi_tsnr_filter(derivRoot, 'RoiDir', D, 'MinTSNR', [45 55])  % absolute
%
% Writes filtered copies of the left/right masks, rebuilds bilat as the union
% of the two filtered sides, and writes the tSNR image it used so you can look
% at where the dropout is.
%
% ---------------------------------------------------------------------------
%  PARAMETERS - the ones you will actually change
% ---------------------------------------------------------------------------
%   'KeepPct'   [left right] percent of each ROI to keep, by tSNR rank.
%               Default [50 50]. A scalar sets both. THIS IS A RANK RULE: n
%               stays the same across subjects and sessions, which is what you
%               want for a group test.
%   'MinTSNR'   [left right] absolute tSNR cut instead. Overrides KeepPct.
%               Lets n vary between subjects, which weights them unequally -
%               use it only if you have a reason.
%   'Sessions'  default {'ses-01'}. 2 and 3 are post-intervention.
%   'Prefix'    functional prefix. Default 'u' -> usub-*, i.e. unwarped and
%               UNSMOOTHED. Use 'ua' for a tree where slice timing ran first.
%               DO NOT point this at su*: see WHY UNSMOOTHED below.
%   'Combine'   how the three runs are combined into one tSNR per voxel:
%               'mean' (default) or 'min'. 'min' is the stricter choice and
%               arguably the right one, because SPM's implicit mask is an
%               INTERSECTION across runs, so the worst run already decides
%               whether a voxel survives at all.
%   'DryRun'    true = report the sweep table and write nothing. Start here.
%
%  The rest, which have sensible defaults:
%   'RoiDir'    flat dir of r<sub>_<ses>_hemi-*_label-*_mask.nii.
%               Default <derivRoot>/../GLM_runauto/ROI
%   'OutDir'    default <RoiDir>_tsnr_L<pct>_R<pct>, so THE DIRECTORY NAME IS
%               THE RECORD of what was done. Filenames are unchanged, so
%               nothing in config.sh needs editing - just point the ROI dir
%               at the new folder.
%   'Subjects'  {} = all
%   'Tasks'     default [1 2 3], matching TASKS in config.sh
%   'HPFSec'    high-pass in seconds, default 128, matching GLM.HPF_cutoff
%   'TR'        default 2.06. Read from the header when it looks sane.
%   'MinVoxels' refuse to write a mask below this, default 5
%   'CSV'       default <OutDir>/roi_tsnr.csv, plus _sweep.csv
%
% ---------------------------------------------------------------------------
%  WHY UNSMOOTHED
% ---------------------------------------------------------------------------
% At FWHM 5 mm on 2 x 2 x 3.45 mm voxels, smoothing multiplies tSNR by 5.4 -
% a variance-equivalent average over ~30 voxels. princeErC is 28-46 voxels per
% hemisphere, so the kernel is as wide as the structure: smoothed tSNR in the
% ErC measures how much tissue surrounds a voxel, which near the collateral
% sulcus is mostly distance to the nearest air. It is not a property of the
% voxel. It also makes neighbouring values nearly equal, which collapses the
% ranking inside a 40-voxel ROI into "pick the middle of the ROI".
%
% So select unsmoothed. If you then analyse smoothed data, the number you
% selected on does not describe the voxel you test - keep both unsmoothed.
%
% ---------------------------------------------------------------------------
%  HOW tSNR IS DEFINED HERE
% ---------------------------------------------------------------------------
%   numerator   mean over time of the RAW timeseries
%   denominator SD of the residual after removing a DCT basis with the GLM's
%               own cutoff - so it is the noise the GLM actually sees, not
%               scanner drift the GLM removes anyway. Same choice as
%               tools/tsnr_report.sh (which uses an FSL Gaussian filter; a DCT
%               set is what SPM itself uses, so this matches the GLM closer).
%
% tSNR is DIRECTION-BLIND, which is the whole point: selecting on it cannot
% bias the GLM2 contrast. Two things it is not blind to, both small and
% pointing opposite ways - the SD includes task-driven variance, so a strongly
% responding voxel reads slightly LOWER; and if you used ResMS instead, the
% grid regressors are in that model, so a responding voxel would read slightly
% HIGHER. Never select on voxelwiseAmplitude: that one is direction-dependent
% and is what GLM2 tests.
%
% ---------------------------------------------------------------------------
%  PICK THE THRESHOLD ON REPRODUCIBILITY, NOT ON THE CONTRAST
% ---------------------------------------------------------------------------
% Run with 'DryRun' true, read the sweep table, and choose. Then judge the
% choice by whether phi becomes more repeatable - run-to-run spread, or median
% R-runs from tools/grid_maps - NOT by whether the contrast grew. The contrast
% is cos(6*dphi) with a sign flip at 15 degrees, so sweeping thresholds and
% keeping the best one finds significance in noise reliably.
%
% ONE SUBJECT x SESSION PER CALL IS THE UNIT OF PARALLELISM - see
% run_roi_tsnr_filter.sh.
%
% THIS WHOLE DIRECTORY IS OPTIONAL. Nothing in the pipeline calls it.

p = inputParser;
p.addParameter('RoiDir',    '');
p.addParameter('OutDir',    '');
p.addParameter('Subjects',  {});
p.addParameter('Sessions',  {'ses-01'});
p.addParameter('KeepPct',   [50 50]);
p.addParameter('MinTSNR',   []);
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

keepPct = o.KeepPct;
if isscalar(keepPct), keepPct = [keepPct keepPct]; end
minTsnr = o.MinTSNR;
if isscalar(minTsnr), minTsnr = [minTsnr minTsnr]; end
useAbs = ~isempty(minTsnr);
if any(keepPct <= 0 | keepPct > 100)
    error('KeepPct must be in (0, 100]; got [%g %g]', keepPct(1), keepPct(2));
end

roiDir = o.RoiDir;
if isempty(roiDir)
    roiDir = fullfile(fileparts(derivRoot), 'GLM_runauto', 'ROI');
end
if ~isfolder(roiDir), error('No such ROI directory: %s', roiDir); end

outDir = o.OutDir;
if isempty(outDir)
    if useAbs
        outDir = sprintf('%s_tsnr_Labs%g_Rabs%g', roiDir, minTsnr(1), minTsnr(2));
    else
        outDir = sprintf('%s_tsnr_L%g_R%g', roiDir, keepPct(1), keepPct(2));
    end
end
if ~o.DryRun && ~isfolder(outDir), mkdir(outDir); end
csvPath = o.CSV;
if isempty(csvPath), csvPath = fullfile(outDir, 'roi_tsnr.csv'); end
sweepPath = regexprep(csvPath, '\.csv$', '_sweep.csv');

SWEEP_PCTS = [100 75 50 25];   % the dry-run table

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
fprintf('  Functionals : %s  (%ssub-*, UNSMOOTHED)\n', derivRoot, o.Prefix);
fprintf('  ROI masks   : %s\n', roiDir);
if useAbs
    fprintf('  Rule        : tSNR >= %g (left), >= %g (right)   ABSOLUTE\n', minTsnr);
else
    fprintf('  Rule        : keep top %g%% (left), top %g%% (right)   BY RANK\n', keepPct);
end
fprintf('  Runs        : tasks [%s], combined by %s\n', num2str(o.Tasks), lower(o.Combine));
fprintf('  High-pass   : DCT, %g s\n', o.HPFSec);
if o.DryRun
    fprintf('  DRY RUN - reporting only, nothing will be written.\n');
else
    fprintf('  Out         : %s\n', outDir);
end
fprintf('============================================\n\n');

rows = {}; sweepRows = {};
for k = 1:size(keys, 1)
    sub = keys{k,1}; ses = keys{k,2}; label = keys{k,3};

    funcs = find_runs(derivRoot, sub, ses, o.Prefix, o.Tasks);
    if isempty(funcs)
        fprintf('%-11s %-8s  no %ssub-* runs found - skipped\n', sub, ses, o.Prefix);
        continue;
    end

    % The 3D mask header is the write reference: same grid as the functionals
    % (checked below), and niftiwrite wants a 3D info for a 3D volume.
    refPath = fullfile(roiDir, sprintf('r%s_%s_hemi-left_label-%s_mask.nii', sub, ses, label));
    info    = niftiinfo(refPath);

    [tsnr, funcInfo, brainTsnr] = tsnr_volume(funcs, o.HPFSec, o.TR, lower(o.Combine));
    check_same_grid(info, funcInfo, refPath, funcs{1});

    fprintf('%s %s  label-%s  %d runs  brain tSNR %.1f\n', ...
            sub, ses, label, numel(funcs), brainTsnr);

    sides = {'left', 'right'};
    kept  = struct('left', [], 'right', []);
    for h = 1:2
        side  = sides{h};
        mPath = fullfile(roiDir, sprintf('r%s_%s_hemi-%s_label-%s_mask.nii', ...
                                         sub, ses, side, label));
        if ~isfile(mPath)
            fprintf('   %-5s  no mask (%s) - skipped\n', side, mPath);
            continue;
        end
        M = niftiread(mPath);
        if numel(M) ~= numel(tsnr)
            error(['ROI mask and the functionals are on different grids:\n' ...
                   '  %s (%d voxels)\n  functional (%d voxels)'], ...
                   mPath, numel(M), numel(tsnr));
        end
        roiIdx = find(M(:) > 0.5);
        tv     = tsnr(roiIdx);
        ok     = isfinite(tv);

        % the sweep table: what each threshold would keep
        fprintf('   %-5s  %3d voxels, tSNR median %.1f (range %.1f-%.1f)\n', ...
                side, numel(roiIdx), median(tv(ok)), min(tv(ok)), max(tv(ok)));
        for q = SWEEP_PCTS
            [sel, cut] = pick_by_pct(tv, ok, q);
            fprintf('            top %3g%% -> %3d voxels, tSNR >= %.1f\n', ...
                    q, sum(sel), cut);
            sweepRows = [sweepRows; {o.Trial, sub, ses, side, label, q, ...
                                     numel(roiIdx), sum(sel), cut}];  %#ok<AGROW>
        end

        if useAbs
            sel = ok & tv >= minTsnr(h);
            cut = minTsnr(h);
        else
            [sel, cut] = pick_by_pct(tv, ok, keepPct(h));
        end

        rows = [rows; {o.Trial, sub, ses, side, label, numel(funcs), ...
                       numel(roiIdx), sum(sel), ...
                       100 * sum(sel) / max(numel(roiIdx), 1), cut, ...
                       median(tv(sel), 'omitnan'), ...
                       median(tv(ok & ~sel), 'omitnan'), ...
                       median(tv(ok), 'omitnan'), brainTsnr}];        %#ok<AGROW>

        if sum(sel) < o.MinVoxels
            warning('roi_tsnr_filter:tooFew', ...
                ['%s %s hemi-%s would keep only %d voxels (MinVoxels %d). ' ...
                 'NOT WRITTEN. Loosen the rule, or fix mthresh and use the ' ...
                 'whole ErC before filtering.'], ...
                 sub, ses, side, sum(sel), o.MinVoxels);
            continue;
        end
        kept.(side) = roiIdx(sel);

        if ~o.DryRun
            write_mask(info, size(tsnr), roiIdx(sel), ...
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
    append_csv(sweepPath, sweepRows, ...
        'trial,subject,session,hemi,label,keep_pct,n_roi,n_kept,cut_tsnr', ...
        '%s,%s,%s,%s,%s,%g,%d,%d,%.3f\n');
    write_info(outDir, derivRoot, roiDir, o, keepPct, minTsnr, useAbs);
    fprintf('\nmasks + tSNR images -> %s\n', outDir);
    fprintf('%d rows -> %s\n', size(rows, 1), csvPath);
    fprintf('%d rows -> %s\n', size(sweepRows, 1), sweepPath);
    fprintf('\nTo use them: point the ROI directory at the folder above.\n');
    fprintf('Filenames are unchanged, so config.sh needs no edit.\n');
else
    fprintf('\nDry run - nothing written. Pick a rule from the table above.\n');
end
fprintf(['Judge the choice by whether phi gets more repeatable ' ...
         '(grid_maps R-runs), not by the contrast.\n']);

end

% ======================================================================

function funcs = find_runs(derivRoot, sub, ses, prefix, tasks)
% One 4D file per task. The glob is anchored on the prefix so that 'u' cannot
% also pick up 'ua' or 'su' files, which coexist in the same folder.
funcs = {};
fdir = fullfile(derivRoot, sub, ses, 'func');
for t = tasks(:)'
    g = dir(fullfile(fdir, sprintf('%ssub-*_%s_task-%d*_bold.nii', prefix, ses, t)));
    g = g(~startsWith({g.name}, 'mean'));
    for a = 1:numel(g)
        nm = g(a).name;
        % reject a longer prefix: the char before "sub-" must end the prefix
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
    brain = mu > mean(mu) / 8;                      % same rule as tsnr_report.sh

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
% divide by the runs ACTUALLY used, which is not numel(funcs) if one was skipped
if strcmp(how, 'mean'), acc = acc / nUsed; end
tsnr = acc;
funcInfo = info;
brainTsnr = median(tsnr(isfinite(tsnr)));
end

function check_same_grid(maskInfo, funcInfo, maskPath, funcPath)
% The resliced masks are written into the EPI grid, so these must agree. If
% they do not, the ROI is being read at the wrong place and every number below
% would be quietly wrong.
a = maskInfo.Transform.T; b = funcInfo.Transform.T;
if ~isequal(maskInfo.ImageSize(1:3), funcInfo.ImageSize(1:3)) || ~all(abs(a(:) - b(:)) < 1e-3)
    error(['Mask and functional are not on the same grid:\n  %s\n  %s\n' ...
           'Use the resliced r*_mask.nii from the same derivatives tree as ' ...
           'the functionals.'], maskPath, funcPath);
end
end

function C = dct_basis(n, k)
% Orthonormal DCT-II basis, column 1 constant - the same set SPM's own
% high-pass removes. Orthonormal, so the residual is Y - C*(C'*Y).
k = max(min(k, n), 1);
C = zeros(n, k);
C(:, 1) = 1 / sqrt(n);
j = 2:k;
C(:, j) = sqrt(2 / n) * cos(pi * (2 * (0:n-1)' + 1) * (j - 1) / (2 * n));
end

function [sel, cut] = pick_by_pct(tv, ok, pct)
% Keep the top pct% by tSNR. Ties are broken by taking them all, so n can
% come out slightly above the nominal count - better than an arbitrary split.
sel = false(size(tv));
nOk = sum(ok);
if nOk == 0, cut = NaN; return; end
nKeep = max(1, round(nOk * pct / 100));
v = sort(tv(ok), 'descend');
cut = v(min(nKeep, numel(v)));
sel = ok & tv >= cut;
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
info.raw.datatype          = 16;      % keep raw in step with Datatype
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

function write_info(outDir, derivRoot, roiDir, o, keepPct, minTsnr, useAbs)
% The directory name carries the rule; this carries everything else, so a
% folder of masks is never anonymous.
fid = fopen(fullfile(outDir, 'filter_info.txt'), 'w');
if fid < 0, return; end
fprintf(fid, 'roi_tsnr_filter\n');
fprintf(fid, 'written          %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, 'functionals      %s (prefix %ssub-*, unsmoothed)\n', derivRoot, o.Prefix);
fprintf(fid, 'source ROI dir   %s\n', roiDir);
if useAbs
    fprintf(fid, 'rule             tSNR >= %g (left), >= %g (right)\n', minTsnr);
else
    fprintf(fid, 'rule             top %g%% (left), top %g%% (right) by tSNR rank\n', keepPct);
end
fprintf(fid, 'tasks            [%s]\n', num2str(o.Tasks));
fprintf(fid, 'runs combined by %s\n', lower(o.Combine));
fprintf(fid, 'high-pass        DCT, %g s\n', o.HPFSec);
fprintf(fid, 'min voxels       %d\n', o.MinVoxels);
fprintf(fid, 'trial            %s\n', o.Trial);
fprintf(fid, '\nbilat masks are the union of the FILTERED left and right.\n');
fclose(fid);
end
