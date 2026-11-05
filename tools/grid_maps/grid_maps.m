function grid_maps(glmRoot, varargin)
%GRID_MAPS  Voxel-wise maps from a GridCAT output directory. No new GLM.
%
%   grid_maps(glmRoot)
%   grid_maps(glmRoot, 'Subjects', {'sub-01s17'}, 'MaskDir', '/path/ROI')
%   grid_maps(glmRoot, 'Subjects', {'sub-01s17'}, 'AnatRoot', '/path/derivatives/spm-preproc')
%
% GridCAT's ROI enters GLM2 only as the scalar phi, so every image it writes
% is whole-brain. This turns those images into maps that say WHERE to look.
% It reads only what GLM1 already wrote and fits no model.
%
% ONE SUBJECT PER CALL IS THE UNIT OF PARALLELISM - see run_grid_maps.sh.
%
% WHAT IT WRITES, per subject x session, into <glmRoot>/grid_maps/
% ----------------------------------------------------------------
%   R-runs      Does this voxel keep its orientation across runs? The
%               resultant length of the per-run unit vectors, 0-1. 1 = every
%               run points the same way in x*alpha space, 0 = they cancel.
%               NaN where fewer than 2 runs have an estimate.
%   R-runs-w    The same, amplitude-weighted - equals
%               nRuns*A_allRunsAvg/sum(A_run). Comparable with GridCAT's own
%               numbers, but one loud run can carry it, so read R-runs first.
%   n-runs      How many runs contributed. R from 2 runs is not comparable
%               with R from 3, so threshold on this before reading R.
%   coh-local   Local between-voxel coherence: resultant length of the unit
%               orientation vectors inside a sphere of RadiusMM, computed on
%               the across-run mean vector. The per-voxel form of the
%               Rayleigh r in your grid_metrics files.
%   z-local     n-local * coh-local^2 - Rayleigh's statistic, so ~1 is
%               chance and >3 is roughly p<0.05 for one neighbourhood.
%               THIS IS THE "where does Rayleigh z cluster" MAP.
%   n-local     Valid voxels in each sphere. An edge voxel with n=6 and a
%               deep one with n=60 are not comparable; this is how you tell.
%   sens        1/sqrt(ResMS). The GLM's own sensitivity map, and the only
%               map here that is BLIND TO DIRECTION - so the only one you
%               may select voxels on without biasing GLM2.
%   roi-kept_*  One per ROI mask: 0 outside the ROI, 1 ROI voxel the GLM
%               dropped, 2 ROI voxel with an estimate. Overlay it to see
%               which part of the ROI the implicit mask deleted and whether
%               it is a contiguous slab (susceptibility) or scattered
%               (noise).
%
% and per subject, when AnatRoot is given and the subject has 2+ sessions:
%
%   R-sessions  Does this voxel keep its orientation ACROSS SESSIONS, in the
%               first session's grid. The environment is identical across
%               your sessions, so this is the strongest within-subject test
%               here - and it does not touch the group design.
%   n-sessions  Sessions contributing at each voxel.
%
% mask.nii IN A GLM DIRECTORY IS NOT AN ROI. It is SPM's implicit analysis
% mask (value > mthresh * global in every scan of every run). This function
% never reads it and never treats it as an ROI. ROIs come only from MaskDir,
% as r<sub>_<ses>_hemi-*_mask.nii, the same convention as the rest of the
% pipeline. "Dropped" above means the GLM has no estimate there, which is
% what the implicit mask having excluded the voxel looks like from here.
%
% ANGLES ARE NEVER AVERAGED. Everything works on the components
% (beta_cos, beta_sin), recovered exactly from the two images GridCAT writes:
%   beta_cos = A*cos(x*ori*pi/180),  beta_sin = A*sin(x*ori*pi/180)
% because ori = atan2(beta_sin, beta_cos)/x and A = hypot(beta_sin, beta_cos).
% Verified to reproduce GridCAT's own allRunsAvg images to float32 precision.
% Summing, filtering and interpolating are all correct on components and all
% wrong on the degree image.
%
% THE NULL FOR R IS NOT ZERO. For n vectors with uniform random phases E[R]
% is about 0.71 at n=2 and 0.56 at n=3 - the reference level this prints.
% And because all runs replay the SAME ROUTE, a voxel can score high from
% route-locked artefact rather than grid coding. Only a permutation null
% (shuffle direction labels within run, refit GLM1, recompute) separates
% those, and that is outside this tool.
%
% IF THE FUNCTIONALS WERE SMOOTHED, z-local IS NOT INTERPRETABLE. A kernel
% makes neighbouring voxels share signal by construction, so the map lights
% up wherever there is any signal at all. Run it on an unsmoothed trial.
% Running it at both and asking whether the smoothed map is merely a blurred
% copy is itself the cleanest test of that confound.
%
% OPTIONS
%   'Subjects'  {} = all. One per call for parallelism.
%   'Sessions'  {} = all
%   'MaskDir'   ROI masks. Empty = skip roi-kept and the CSV.
%   'AnatRoot'  preprocessed derivatives tree holding
%               sub-*/ses-*/anat/sub-*_T2w.nii. Empty (default) = skip the
%               cross-session map, and then SPM is not needed at all.
%               MUST be the tree the GLM output was built from: each
%               session's T2w sits in that session's EPI world space, and
%               this checks the resliced ROI mask's affine to catch a
%               mismatched tree.
%   'OutDir'    default <glmRoot>/grid_maps
%   'Event'     grid event name, default 'translate'
%   'XFold'     default 6
%   'RadiusMM'  sphere radius for coh-local in MILLIMETRES, default 5.
%               Built from the voxel sizes, so anisotropy is handled - at
%               2x2x3.45 mm that is +/-2 voxels in x,y and +/-1 in z.
%   'Trial'     label written into the CSV so trials stack up
%   'CSV'       default <OutDir>/grid_maps.csv
%
% THIS WHOLE DIRECTORY IS OPTIONAL. Nothing in the pipeline calls it.

p = inputParser;
p.addParameter('Subjects', {});
p.addParameter('Sessions', {});
p.addParameter('MaskDir',  '');
p.addParameter('AnatRoot', '');
p.addParameter('OutDir',   '');
p.addParameter('Event',    'translate');
p.addParameter('XFold',    6);
p.addParameter('RadiusMM', 5);
p.addParameter('Trial',    '');
p.addParameter('CSV',      '');
p.parse(varargin{:});
o = p.Results;

if ~isfolder(glmRoot), error('No such GLM output directory: %s', glmRoot); end
xf = o.XFold;

outDir = o.OutDir;
if isempty(outDir), outDir = fullfile(glmRoot, 'grid_maps'); end
if ~isfolder(outDir), mkdir(outDir); end
csvPath = o.CSV;
if isempty(csvPath), csvPath = fullfile(outDir, 'grid_maps.csv'); end

if ~isempty(o.AnatRoot)
    if ~isfolder(o.AnatRoot), error('No such AnatRoot: %s', o.AnatRoot); end
    if isempty(which('spm_coreg'))
        error(['AnatRoot was given but SPM is not on the MATLAB path. ' ...
               'The cross-session map needs spm_coreg to register the ' ...
               'sessions - they are NOT in a common space (one T2w per ' ...
               'session, each coregistered only to its own EPI).']);
    end
end

% ---- discover subject x session ---------------------------------------
glm1 = dir(fullfile(glmRoot, '*_GLM1'));
glm1 = glm1([glm1.isdir]);
subs = {}; sess = {}; dirs = {};
for a = 1:numel(glm1)
    tok = regexp(glm1(a).name, '^(sub-[^_]+)_(ses-[^_]+)_GLM1$', 'tokens', 'once');
    if isempty(tok), continue; end
    if ~isempty(o.Subjects) && ~ismember(tok{1}, o.Subjects), continue; end
    if ~isempty(o.Sessions) && ~ismember(tok{2}, o.Sessions), continue; end
    subs{end+1} = tok{1}; sess{end+1} = tok{2};                   %#ok<AGROW>
    dirs{end+1} = fullfile(glm1(a).folder, glm1(a).name);         %#ok<AGROW>
end
if isempty(subs), error('No matching *_GLM1 directories under %s', glmRoot); end
[subs, ord] = sort(subs); sess = sess(ord); dirs = dirs(ord);
uSubs = unique(subs, 'stable');

fprintf('\n============================================\n');
fprintf('  Grid maps\n');
fprintf('  GLM1 output : %s\n', glmRoot);
fprintf('  ROI masks   : %s\n', pick(isempty(o.MaskDir),  '<none - skipping roi-kept>', o.MaskDir));
fprintf('  Anat tree   : %s\n', pick(isempty(o.AnatRoot), '<none - skipping R-sessions>', o.AnatRoot));
fprintf('  Maps out    : %s\n', outDir);
fprintf('  %d-fold, local sphere radius %g mm\n', xf, o.RadiusMM);
fprintf('  R under uniform phases: E[R] = %.2f (n=2), %.2f (n=3)\n', uniformR(2), uniformR(3));
fprintf('  z-local is NOT interpretable on smoothed functionals.\n');
fprintf('============================================\n\n');
fprintf('%-11s %-8s %5s %-26s %6s %6s %6s %7s %7s\n', ...
        'subject','session','n','roi','n_roi','n_kept','%kept','medR','med_z');

rows = {};
for s = 1:numel(uSubs)
    sub  = uSubs{s};
    idx  = find(strcmp(subs, sub));
    keep = struct('ses', {}, 'mc', {}, 'ms', {}, 'info', {}, 'dir', {});

    for b = idx
        ses  = sess{b};
        gdir = dirs{b};
        [bc, bs, runNums, refInfo] = load_components(gdir, o.Event, xf);
        if isempty(runNums)
            fprintf('%-11s %-8s  (no %s run images - skipped)\n', sub, ses, o.Event);
            continue;
        end
        pd = double(refInfo.PixelDimensions(1:3));

        % ---- across-run agreement, per voxel --------------------------
        A     = hypot(bc, bs);
        valid = isfinite(bc) & isfinite(bs) & A > 0;   % A==0 is the atan2(0,0) trap
        nRuns = sum(valid, 4);

        % unweighted: unit vectors, so amplitude cannot carry the result
        [Su, Sv] = sum_valid(bc ./ A, bs ./ A, valid);
        Rruns = hypot(Su, Sv) ./ max(nRuns, 1);
        Rruns(nRuns < 2) = NaN;

        % weighted, and the across-run mean vector (== GridCAT's allRunsAvg)
        [Sbc, Sbs] = sum_valid(bc, bs, valid);
        SA = sum_valid(A, [], valid);
        RrunsW = hypot(Sbc, Sbs) ./ SA;
        RrunsW(nRuns < 2 | ~(SA > 0)) = NaN;

        mc = Sbc ./ max(nRuns, 1);  mc(nRuns < 1) = NaN;
        ms = Sbs ./ max(nRuns, 1);  ms(nRuns < 1) = NaN;

        % ---- local between-voxel coherence, on the mean vector --------
        Am = hypot(mc, ms);
        vm = isfinite(Am) & Am > 0;
        [Lu, N] = sphere_sum(mc ./ Am, vm, pd, o.RadiusMM);
        Lv      = sphere_sum(ms ./ Am, vm, pd, o.RadiusMM);
        coh  = hypot(Lu, Lv) ./ max(N, 1);
        zloc = N .* coh.^2;
        thin = N < 5;                       % too few members for a Rayleigh
        coh(thin) = NaN;  zloc(thin) = NaN;

        % ---- direction-blind sensitivity ------------------------------
        sens    = [];
        resPath = fullfile(gdir, 'ResMS.nii');
        if isfile(resPath)
            res  = double(niftiread(resPath));
            sens = 1 ./ sqrt(res);
            sens(~isfinite(res) | res <= 0) = NaN;
        end

        stem = fullfile(outDir, sprintf('%s_%s', sub, ses));
        write_map(refInfo, Rruns,  [stem '_R-runs']);
        write_map(refInfo, RrunsW, [stem '_R-runs-w']);
        write_map(refInfo, nRuns,  [stem '_n-runs']);
        write_map(refInfo, coh,    [stem '_coh-local']);
        write_map(refInfo, zloc,   [stem '_z-local']);
        write_map(refInfo, N,      [stem '_n-local']);
        if ~isempty(sens), write_map(refInfo, sens, [stem '_sens']); end

        keep(end+1) = struct('ses', ses, 'mc', mc, 'ms', ms, ...
                             'info', refInfo, 'dir', gdir);        %#ok<AGROW>

        maps = struct('R_runs', Rruns, 'coh', coh, 'z', zloc, ...
                      'sens', sens, 'amp', Am, 'kept', vm);
        rows = [rows; roi_rows(o, sub, ses, numel(runNums), maps, ...
                               refInfo, stem, 'roi-kept')];       %#ok<AGROW>
    end

    % ---- cross-session map --------------------------------------------
    if isempty(o.AnatRoot) || numel(keep) < 2, continue; end
    [Rses, nSes, tgtInfo] = sessions_map(o, sub, keep);
    if isempty(Rses), continue; end
    stem = fullfile(outDir, sprintf('%s_across-sessions', sub));
    write_map(tgtInfo, Rses, [stem '_R-sessions']);
    write_map(tgtInfo, nSes, [stem '_n-sessions']);

    maps = struct('R_runs', Rses, 'coh', nan(size(Rses)), 'z', nan(size(Rses)), ...
                  'sens', [], 'amp', nan(size(Rses)), 'kept', nSes >= 2);
    rows = [rows; roi_rows(o, sub, keep(1).ses, numel(keep), maps, ...
                           tgtInfo, stem, 'roi-across')];          %#ok<AGROW>
end

if ~isempty(rows)
    append_csv(csvPath, rows);
    fprintf('\n%d rows -> %s\n', size(rows, 1), csvPath);
end
fprintf('maps -> %s\n', outDir);

end

% ======================================================================
%  loading
% ======================================================================

function [bc, bs, runNums, refInfo] = load_components(gdir, event, xf)
% Recover (beta_cos, beta_sin) exactly from the two images GridCAT writes.
f = dir(fullfile(gdir, sprintf('voxelwiseGridOri_%s_run*_deg.nii', event)));
runNums = []; bc = []; bs = []; refInfo = [];
for k = 1:numel(f)
    n = regexp(f(k).name, '_run(\d+)_deg\.nii$', 'tokens', 'once');
    if isempty(n), continue; end
    r     = str2double(n{1});
    aPath = fullfile(gdir, sprintf('voxelwiseAmplitude_%s_run%d.nii', event, r));
    if ~isfile(aPath)
        warning('grid_maps:noAmplitude', 'No amplitude image for run %d in %s', r, gdir);
        continue;
    end
    if isempty(refInfo), refInfo = niftiinfo(fullfile(gdir, f(k).name)); end
    ori = double(niftiread(fullfile(gdir, f(k).name)));   % degrees, [0, 360/xf)
    A   = double(niftiread(aPath));
    ang = xf * ori * pi / 180;                            % back to the raw angle
    % explicit index: size([],4) is 1, so end+1 would leave slice 1 as zeros
    nk = numel(runNums) + 1;
    runNums(nk)  = r;                        %#ok<AGROW>
    bc(:,:,:,nk) = A .* cos(ang);            %#ok<AGROW>
    bs(:,:,:,nk) = A .* sin(ang);            %#ok<AGROW>
end
[runNums, ord] = sort(runNums);
if ~isempty(runNums), bc = bc(:,:,:,ord); bs = bs(:,:,:,ord); end
end

function [Sa, Sb] = sum_valid(a, b, valid)
% Sum over the 4th dimension counting only valid entries.
a(~valid) = 0;  Sa = sum(a, 4);
Sb = [];
if ~isempty(b), b(~valid) = 0; Sb = sum(b, 4); end
end

% ======================================================================
%  local coherence
% ======================================================================

function [S, N] = sphere_sum(X, valid, pd, radiusMM)
% Sum X and count valid voxels over a sphere of radiusMM, zero-padded at the
% edges (never wrapped). The radius is in mm and the index range comes from
% the voxel sizes, so anisotropic voxels give an anisotropic range.
X(~valid) = 0;
S = zeros(size(X));
N = zeros(size(X));
V = double(valid);
pd = pd(:)';
rv = floor(radiusMM ./ pd);
for dx = -rv(1):rv(1)
    for dy = -rv(2):rv(2)
        for dz = -rv(3):rv(3)
            if sum(([dx dy dz] .* pd).^2) > radiusMM^2, continue; end
            S = S + shift3(X, dx, dy, dz);
            N = N + shift3(V, dx, dy, dz);
        end
    end
end
end

function Y = shift3(X, dx, dy, dz)
% Y(i,j,k) = X(i+dx, j+dy, k+dz), zeros outside. Not circshift: wrapping
% would join opposite edges of the field of view.
[nx, ny, nz] = size(X);
Y = zeros(nx, ny, nz);
i1 = max(1, 1-dx); i2 = min(nx, nx-dx);
j1 = max(1, 1-dy); j2 = min(ny, ny-dy);
k1 = max(1, 1-dz); k2 = min(nz, nz-dz);
if i1 > i2 || j1 > j2 || k1 > k2, return; end
Y(i1:i2, j1:j2, k1:k2) = X(i1+dx:i2+dx, j1+dy:j2+dy, k1+dz:k2+dz);
end

% ======================================================================
%  cross-session
% ======================================================================

function [Rses, nSes, tgtInfo] = sessions_map(o, sub, keep)
% Resultant length of the per-session mean vectors, in session 1's grid.
%
% The sessions are NOT in a common space: the pipeline gives each session its
% own T2w and coregisters it to that session's own reference EPI, so each
% session's images live in that session's scanner space. Reslicing by header
% alone would silently compare different tissue. So the between-session rigid
% transform is measured here, T2w to T2w, with spm_coreg.
Rses = []; nSes = []; tgtInfo = keep(1).info;
tgtM    = sform_of(tgtInfo);
tgtSize = double(tgtInfo.ImageSize(1:3));

check_anat_matches(o, sub, keep(1).ses, tgtM);

uc = nan([tgtSize numel(keep)]);
us = nan([tgtSize numel(keep)]);
for k = 1:numel(keep)
    if k == 1
        mc = keep(k).mc;  ms = keep(k).ms;      % already the target grid
    else
        W = session_transform(o.AnatRoot, sub, keep(1).ses, keep(k).ses);
        if isempty(W), continue; end
        T  = (W * sform_of(keep(k).info)) \ tgtM;   % target voxel -> source voxel
        mc = reslice_lin(keep(k).mc, T, tgtSize);
        ms = reslice_lin(keep(k).ms, T, tgtSize);
    end
    Am = hypot(mc, ms);
    ok = isfinite(Am) & Am > 0;
    c = mc ./ Am; c(~ok) = NaN;  uc(:,:,:,k) = c;
    d = ms ./ Am; d(~ok) = NaN;  us(:,:,:,k) = d;
end

valid = isfinite(uc) & isfinite(us);
nSes  = sum(valid, 4);
[Su, Sv] = sum_valid(uc, us, valid);
Rses = hypot(Su, Sv) ./ max(nSes, 1);
Rses(nSes < 2) = NaN;
fprintf('%-11s %-8s %5d %-26s\n', sub, 'across', numel(keep), ...
        sprintf('R-sessions, %d vox with 2+', sum(nSes(:) >= 2)));
end

function W = session_transform(anatRoot, sub, sesTgt, sesSrc)
% World-to-world rigid transform taking sesSrc's space into sesTgt's.
% spm_run_coreg writes the source's new mat as spm_matrix(x)\M, so the point
% at world w moves to inv(spm_matrix(x))*w.
W = [];
tgt = find_t2w(anatRoot, sub, sesTgt);
src = find_t2w(anatRoot, sub, sesSrc);
if isempty(tgt) || isempty(src)
    warning('grid_maps:noT2w', ...
            'No T2w for %s %s or %s under %s - skipping that session.', ...
            sub, sesTgt, sesSrc, anatRoot);
    return;
end
fprintf('   coreg %s -> %s ...\n', sesSrc, sesTgt);
x = spm_coreg(spm_vol(tgt), spm_vol(src));
W = inv(spm_matrix(x(:)'));
end

function f = find_t2w(anatRoot, sub, ses)
% The session's own T2w, after preprocessing coregistered it to that
% session's reference EPI - so its header is in that session's world space.
% Anchored on sub- so the r-prefixed resliced masks cannot match.
d = dir(fullfile(anatRoot, sub, ses, 'anat', sprintf('%s_%s_*T2w.nii', sub, ses)));
f = '';
if ~isempty(d), f = fullfile(d(1).folder, d(1).name); end
end

function check_anat_matches(o, sub, ses, glmM)
% The resliced ROI mask in the anat tree sits on the session's reference EPI
% grid. If that does not match the GLM output's grid, the anat tree is from a
% different preprocessing variant and its T2w headers are in the wrong space.
d = dir(fullfile(o.AnatRoot, sub, ses, 'anat', sprintf('r%s_%s_hemi-*_mask.nii', sub, ses)));
if isempty(d), return; end
m = sform_of(niftiinfo(fullfile(d(1).folder, d(1).name)));
if max(abs(m(:) - glmM(:))) > 1e-3
    error(['AnatRoot does not match this GLM output.\n' ...
           '  %s\n is on a different grid from the GLM1 images for %s %s.\n' ...
           'Point AnatRoot at the derivatives tree this GLM output was built from.'], ...
           d(1).name, sub, ses);
end
end

function Y = reslice_lin(X, T, tgtSize)
% Trilinear, on components only. T maps TARGET 0-based voxel coordinates to
% SOURCE 0-based voxel coordinates. NaN corners drop out of the weighted mean
% and the remaining weight is renormalised, so a voxel needs more than half
% its weight from real data - otherwise it stays NaN rather than being
% invented from one neighbour.
[i, j, k] = ndgrid(0:tgtSize(1)-1, 0:tgtSize(2)-1, 0:tgtSize(3)-1);
n = numel(i);
c = T * [i(:)'; j(:)'; k(:)'; ones(1, n)];
x = c(1,:)'; y = c(2,:)'; z = c(3,:)';
[nx, ny, nz] = size(X);
x0 = floor(x); y0 = floor(y); z0 = floor(z);
fx = x - x0;  fy = y - y0;  fz = z - z0;
num = zeros(n, 1); den = zeros(n, 1);
for a = 0:1
    wa = a * fx + (1 - a) * (1 - fx);
    for b = 0:1
        wb = b * fy + (1 - b) * (1 - fy);
        for c3 = 0:1
            w  = wa .* wb .* (c3 * fz + (1 - c3) * (1 - fz));
            xi = x0 + a; yi = y0 + b; zi = z0 + c3;
            in = xi >= 0 & xi < nx & yi >= 0 & yi < ny & zi >= 0 & zi < nz & w > 0;
            if ~any(in), continue; end
            v  = X(sub2ind([nx ny nz], xi(in) + 1, yi(in) + 1, zi(in) + 1));
            ok = isfinite(v);
            ii = find(in); ii = ii(ok);
            num(ii) = num(ii) + w(ii) .* v(ok);
            den(ii) = den(ii) + w(ii);
        end
    end
end
Y = nan(n, 1);
good = den > 0.5;
Y(good) = num(good) ./ den(good);
Y = reshape(Y, tgtSize);
end

function M = sform_of(info)
% NIfTI sform: 0-based voxel -> world. SPM always writes it.
r = info.raw;
M = double([r.srow_x; r.srow_y; r.srow_z; 0 0 0 1]);
if ~any(M(1:3, 1:3) ~= 0)
    error('This image has no sform; cannot place it in world space.');
end
end

% ======================================================================
%  ROIs and output
% ======================================================================

function rows = roi_rows(o, sub, ses, nUnits, maps, refInfo, stem, tag)
% ROIs come ONLY from MaskDir. The GLM's own mask.nii is SPM's implicit
% analysis mask, not an ROI, and is never read.
rows = {};
if isempty(o.MaskDir)
    fprintf('%-11s %-8s %5d %-26s\n', sub, ses, nUnits, '<no MaskDir>');
    return;
end
masks = dir(fullfile(o.MaskDir, sprintf('r%s_%s_hemi-*_mask.nii', sub, ses)));
if isempty(masks)
    fprintf('%-11s %-8s %5d %-26s\n', sub, ses, nUnits, '(no ROI masks)');
    return;
end
for m = 1:numel(masks)
    M = niftiread(fullfile(masks(m).folder, masks(m).name));
    if numel(M) ~= numel(maps.kept)
        error(['ROI mask and the maps are on different grids:\n' ...
               '  %s (%d voxels)\n  map (%d voxels)'], ...
               masks(m).name, numel(M), numel(maps.kept));
    end
    roi  = M(:) > 0.5;
    kept = roi & maps.kept(:);
    label = regexprep(masks(m).name, '^r.*_(hemi-[^_]+_label-[^_]+)_mask\.nii$', '$1');

    K = zeros(size(maps.kept), 'single');
    K(roi) = 1; K(kept) = 2;
    write_map(refInfo, K, sprintf('%s_%s_%s', stem, tag, label));

    nroi = sum(roi); nkept = sum(kept);
    rows = [rows; {o.Trial, sub, ses, label, nUnits, nroi, nkept, ...
                   100 * nkept / max(nroi, 1), ...
                   median(maps.R_runs(kept), 'omitnan'), ...
                   median(maps.coh(kept),    'omitnan'), ...
                   median(maps.z(kept),      'omitnan'), ...
                   med_or_nan(maps.sens, kept), ...
                   median(maps.amp(kept),    'omitnan')}];          %#ok<AGROW>
    fprintf('%-11s %-8s %5d %-26s %6d %6d %5.1f%% %7.3f %7.2f\n', ...
            sub, ses, nUnits, label, nroi, nkept, ...
            100 * nkept / max(nroi, 1), ...
            median(maps.R_runs(kept), 'omitnan'), ...
            median(maps.z(kept), 'omitnan'));
end
end

function v = med_or_nan(X, idx)
v = NaN;
if ~isempty(X), v = median(X(idx), 'omitnan'); end
end

function write_map(refInfo, X, stemNoExt)
info = refInfo;
info.Datatype              = 'single';
info.BitsPerPixel          = 32;
info.ImageSize             = size(X);
info.Description           = 'grid_maps';
info.MultiplicativeScaling = 1;
info.AdditiveOffset        = 0;
niftiwrite(single(X), stemNoExt, info);
end

function r = uniformR(n)
% E[R] for n unit vectors with uniform random phases - the reference level
% for R, which is emphatically not zero.
s = rng(0);
th = 2 * pi * rand(20000, n);
r  = mean(hypot(sum(cos(th), 2), sum(sin(th), 2)) / n);
rng(s);
end

function append_csv(path, rows)
hdr = ['trial,subject,session,roi,n_units,n_roi,n_kept,pct_kept,' ...
       'med_R,med_coh_local,med_z_local,med_sens,med_amp'];
newFile = ~isfile(path);
fid = fopen(path, 'a');
if fid < 0, error('Cannot write CSV: %s', path); end
if newFile, fprintf(fid, '%s\n', hdr); end
for i = 1:size(rows, 1)
    fprintf(fid, '%s,%s,%s,%s,%d,%d,%d,%.2f,%.4f,%.4f,%.4f,%.4f,%.4f\n', rows{i,:});
end
fclose(fid);
end

function s = pick(c, a, b)
if c, s = a; else, s = b; end
end
