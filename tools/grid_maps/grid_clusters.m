function grid_clusters(glmRoot, varargin)
%GRID_CLUSTERS  Clusters of voxels that agree with their neighbours on grid
%               orientation. No new GLM, no Rayleigh sphere, no radius.
%
%   grid_clusters(glmRoot)
%   grid_clusters(glmRoot, 'Subjects', {'sub-01s17'}, 'MaskDir', '/path/ROI')
%
% WHY THIS AND NOT A RAYLEIGH MAP
% -------------------------------
% A local Rayleigh z needs a sphere radius, and the answer moves when the
% radius moves. This asks a question with no free spatial parameter: join two
% face-adjacent voxels if their grid orientations agree to within AngleTol,
% then look at the connected components. You get objects with sizes, and a
% size can be tested against a null.
%
% 1 deg AND 59 deg ARE 2 deg APART at 6-fold symmetry, and this never does
% modular arithmetic to get that right. Each voxel becomes a unit vector
%   u = (cos(x*ori), sin(x*ori))
% and two voxels agree when dot(u_i, u_j) >= cos(x*AngleTol). The wrap is
% built into the dot product; there is no branch to get wrong.
%
% ONE SUBJECT PER CALL IS THE UNIT OF PARALLELISM - see run_grid_clusters.sh.
%
% WHAT IT WRITES, per subject x session, into <glmRoot>/grid_clusters/
% -------------------------------------------------------------------
%   nbr-dev     THE ONE TO OPEN FIRST. Per voxel, how far its orientation
%               sits from its face neighbours' , in ORIENTATION DEGREES,
%               0-30. Read it directly: 0 = identical to its neighbours,
%               15 = exactly chance, 30 = maximally opposed (the sign-flip
%               line is at 15). Needs no threshold and no cluster to read.
%   clust-label Integer cluster id, 0 outside. Ids run largest first.
%   clust-size  Every voxel carries its own cluster's size in voxels, so you
%               can threshold at the null's 95th percentile in any viewer and
%               see only clusters bigger than chance.
%   clust-ori   Every voxel carries its cluster's mean orientation, degrees
%               in [0, 360/x). Constant inside a cluster, so clusters read as
%               flat patches of colour and you can see whether neighbouring
%               clusters hold different orientations or the same one.
%
% and two CSVs: one row per cluster in grid_clusters.csv, one row per
% subject x session in grid_clusters_summary.csv.
%
% THE NULL
% --------
% Cluster extent means nothing on its own, so the orientation values are
% shuffled among the same voxels Perms times and the clustering redone. The
% spatial support and the value distribution are held fixed; only the
% arrangement changes. Two numbers come out:
%
%   max_null_p95   the 95th percentile of the largest cluster under shuffling.
%                  A cluster bigger than this is bigger than chance.
%   dev_p          one-sided p for the whole field being smoother than
%                  chance, from the mean neighbour agreement.
%
% WHAT A PASS DOES AND DOES NOT MEAN. It means the orientation field is
% spatially structured. It does NOT mean grid coding: any source of spatial
% autocorrelation produces it.
%
%   - ON SMOOTHED FUNCTIONALS THIS TEST IS GUARANTEED TO PASS AND IS
%     MEANINGLESS. A kernel makes neighbours share signal by construction.
%     Run it on an unsmoothed trial.
%   - Even unsmoothed, EPI carries autocorrelation from reconstruction and
%     from realignment interpolation, so a small effect is expected.
%   - Whether aligned beats misaligned is still GLM2's question, not this
%     one. This says WHERE TO LOOK.
%
% WHY AngleTol DEFAULTS TO 5 AND NOT 15
% -------------------------------------
% Two random orientations agree within AngleTol with probability
% AngleTol/(180/x) - at 6-fold, AngleTol/30. Each voxel has 6 face
% neighbours, so the shuffled graph has mean degree 6*AngleTol/30.
%
% Bond percolation on a simple cubic lattice has a threshold at p = 0.2488.
% Above it a single connected component spans the whole volume IN PURE
% NOISE. At AngleTol = 15 (the sign-flip line, the tempting choice) p = 0.5,
% and both the data and the null give one giant cluster the size of the
% brain: max cluster extent stops being a statistic. At the default 5,
% p = 0.167 and the null stays in small pieces, so a real patch stands out.
%
% The tool prints p and warns past 0.22, and warns again if the largest
% cluster actually holds more than a quarter of the valid voxels. Two columns
% tell percolation from signal: frac_largest, the share of valid voxels in the
% biggest cluster, and max_clust_R, its internal resultant length. A real
% patch is small and has max_clust_R near 1; a percolating chain is enormous
% and has max_clust_R near 0, because it wanders through orientation space one
% tolerable step at a time. dev_deg and dev_p are immune to all of this and
% stay readable at any tolerance, which is the other reason to read nbr-dev
% first.
%
% It is run whole-brain on purpose. If the field clusters in visual cortex
% and not in the ErC, the method works and the ErC is just noisy; if it
% clusters nowhere, something upstream is wrong. Both are worth knowing, and
% neither is visible from inside an ROI.
%
% OPTIONS
%   'Subjects'  {} = all. One per call for parallelism.
%   'Sessions'  default {'ses-01'}. Sessions 2 and 3 are post-intervention.
%   'MaskDir'   ROI masks, as r<sub>_<ses>_hemi-*_mask.nii. Empty = skip the
%               ROI overlap columns. Clustering is whole-brain either way.
%   'Run'       0 (default) = the across-run mean vector. 1/2/3 = that run
%               alone, for asking whether the same clusters come back.
%   'AngleTol'  agreement tolerance in ORIENTATION degrees, default 5.
%               Must be in (0, 180/XFold) - at 6-fold that is (0, 30), and 30
%               would join everything to everything.
%               DO NOT RAISE THIS ABOVE ~7 WITHOUT READING THE NOTE ON
%               PERCOLATION BELOW. The default is not a round number by
%               accident.
%   'MinSize'   smallest cluster written to the CSV, default 10 voxels. The
%               images keep every cluster.
%   'Perms'     shuffles for the null, default 500. 0 skips the null, and
%               then the size images have nothing to be read against.
%   'Event'     grid event name, default 'translate'
%   'XFold'     default 6
%   'Seed'      RNG seed for the null, default 0, so a rerun reproduces.
%   'OutDir'    default <glmRoot>/grid_clusters
%   'Trial'     label written into both CSVs so trials stack up
%   'CSV'       default <OutDir>/grid_clusters.csv
%
% Reads only what GLM1 wrote, by image name, so keepUnusedGridEvents and any
% other change to the beta NUMBERING is irrelevant here.
%
% THIS WHOLE DIRECTORY IS OPTIONAL. Nothing in the pipeline calls it.

p = inputParser;
p.addParameter('Subjects', {});
p.addParameter('Sessions', {'ses-01'});
p.addParameter('MaskDir',  '');
p.addParameter('Run',      0);
p.addParameter('AngleTol', 5);
p.addParameter('MinSize',  10);
p.addParameter('Perms',    500);
p.addParameter('Event',    'translate');
p.addParameter('XFold',    6);
p.addParameter('Seed',     0);
p.addParameter('OutDir',   '');
p.addParameter('Trial',    '');
p.addParameter('CSV',      '');
p.parse(varargin{:});
o = p.Results;

if ~isfolder(glmRoot), error('No such GLM output directory: %s', glmRoot); end
xf = o.XFold;
if ~(o.AngleTol > 0 && o.AngleTol < 180/xf)
    error(['AngleTol must be in (0, %g) orientation degrees at %d-fold ' ...
           'symmetry. %g was given; %g would join every voxel to every ' ...
           'other.'], 180/xf, xf, o.AngleTol, 180/xf);
end

outDir = o.OutDir;
if isempty(outDir), outDir = fullfile(glmRoot, 'grid_clusters'); end
if ~isfolder(outDir), mkdir(outDir); end
csvPath = o.CSV;
if isempty(csvPath), csvPath = fullfile(outDir, 'grid_clusters.csv'); end
sumPath = regexprep(csvPath, '\.csv$', '_summary.csv');
if strcmp(sumPath, csvPath), sumPath = [csvPath '_summary.csv']; end

% agreement threshold on the dot product, and the chance rate of passing it
dotThr = cos(xf * o.AngleTol * pi / 180);
pAgree = o.AngleTol / (180 / xf);
if pAgree > 0.22
    warning('grid_clusters:percolation', ...
        ['AngleTol = %g gives a chance agreement rate of %.3f, at or above ' ...
         'the 0.2488 bond-percolation threshold for a 6-connected lattice. ' ...
         'The SHUFFLED data will produce one cluster spanning the whole ' ...
         'volume, so cluster extent will not discriminate. Use 5-7 instead. ' ...
         'dev_deg and dev_p are unaffected.'], o.AngleTol, pAgree);
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
    subs{end+1} = tok{1}; sess{end+1} = tok{2};                    %#ok<AGROW>
    dirs{end+1} = fullfile(glm1(a).folder, glm1(a).name);          %#ok<AGROW>
end
if isempty(subs), error('No matching *_GLM1 directories under %s', glmRoot); end
[subs, ord] = sort(subs); sess = sess(ord); dirs = dirs(ord);

if o.Run == 0, runLab = 'mean'; else, runLab = sprintf('%d', o.Run); end

fprintf('\n============================================\n');
fprintf('  Grid orientation clusters\n');
fprintf('  GLM1 output : %s\n', glmRoot);
fprintf('  ROI masks   : %s\n', pick(isempty(o.MaskDir), '<none - no ROI columns>', o.MaskDir));
fprintf('  Clusters out: %s\n', outDir);
fprintf('  %d-fold, run %s, agreement within %g deg of orientation\n', xf, runLab, o.AngleTol);
fprintf('  (two voxels join when dot(u_i,u_j) >= %.4f)\n', dotThr);
fprintf('  chance agreement rate %.3f (percolation sets in near 0.249)\n', pAgree);
fprintf('  null: %d shuffles, seed %d\n', o.Perms, o.Seed);
fprintf('  NOT interpretable on smoothed functionals.\n');
fprintf('============================================\n\n');
fprintf('%-11s %-8s %7s %7s %8s %8s %9s %7s %6s %6s\n', ...
        'subject','session','n_vox','dev','dev_null','p','max_clust', ...
        'null95','frac','clustR');

cRows = []; sRows = [];
for b = 1:numel(subs)
    sub = subs{b}; ses = sess{b}; gdir = dirs{b};

    [u, v, valid, info, nRunsUsed] = orientation_field(gdir, o.Event, xf, o.Run);
    if isempty(u)
        fprintf('%-11s %-8s  (no %s images for run %s - skipped)\n', sub, ses, o.Event, runLab);
        continue;
    end

    vi = find(valid);
    nv = numel(vi);
    if nv < 2
        fprintf('%-11s %-8s  (only %d valid voxels - skipped)\n', sub, ses, nv);
        continue;
    end
    uv = u(vi); vv = v(vi);

    % ---- adjacency among valid voxels, 6-connectivity -----------------
    [ia, ib] = face_pairs(valid);
    if isempty(ia)
        fprintf('%-11s %-8s  (no adjacent valid voxels - skipped)\n', sub, ses);
        continue;
    end
    d = uv(ia) .* uv(ib) + vv(ia) .* vv(ib);      % cos of the raw angle gap

    % ---- per-voxel neighbour deviation, in orientation degrees --------
    sumd = accumarray([ia; ib], [d; d], [nv 1]);
    cnt  = accumarray([ia; ib], 1,      [nv 1]);
    md   = sumd ./ max(cnt, 1);
    devV = acosd(min(max(md, -1), 1)) / xf;        % 0 .. 180/xf
    devV(cnt == 0) = NaN;

    devImg = nan(size(valid));
    devImg(vi) = devV;

    % ---- clusters ------------------------------------------------------
    [lab, sz] = components(ia(d >= dotThr), ib(d >= dotThr), nv);

    % relabel largest first, so id 1 is the biggest cluster
    [sz, newOrder] = sort(sz, 'descend');
    remap = zeros(numel(newOrder), 1);
    remap(newOrder) = 1:numel(newOrder);
    lab = remap(lab);

    labImg  = zeros(size(valid), 'single');
    sizeImg = zeros(size(valid), 'single');
    oriImg  = nan(size(valid));
    labImg(vi)  = lab;
    sizeImg(vi) = sz(lab);

    % per-cluster circular mean orientation, amplitude and centroid
    nC   = numel(sz);
    Su   = accumarray(lab, uv, [nC 1]);
    Sv   = accumarray(lab, vv, [nC 1]);
    ang  = atan2(Sv, Su);
    ang(ang < 0) = ang(ang < 0) + 2*pi;
    oriC = ang / xf * 180/pi;                      % [0, 360/xf)
    Rin  = hypot(Su, Sv) ./ sz(:);
    oriImg(vi) = oriC(lab);

    % ---- the null ------------------------------------------------------
    devObs = acosd(min(max(mean(d), -1), 1)) / xf;
    maxNull = []; devNull = []; p95 = NaN; devP = NaN; devNullMean = NaN;
    if o.Perms > 0
        rs = rng(o.Seed);
        maxNull = zeros(o.Perms, 1);
        devNull = zeros(o.Perms, 1);
        for t = 1:o.Perms
            pr = randperm(nv);
            up = uv(pr); vp = vv(pr);
            dp = up(ia) .* up(ib) + vp(ia) .* vp(ib);
            [~, szp] = components(ia(dp >= dotThr), ib(dp >= dotThr), nv);
            maxNull(t) = max(szp);
            devNull(t) = acosd(min(max(mean(dp), -1), 1)) / xf;
        end
        rng(rs);
        p95         = prctile_(maxNull, 95);
        devNullMean = mean(devNull);
        % one-sided: is the field SMOOTHER (smaller deviation) than chance
        devP = (1 + sum(devNull <= devObs)) / (o.Perms + 1);
    end

    % ---- ROI overlap ---------------------------------------------------
    roiMask = []; roiName = '';
    if ~isempty(o.MaskDir)
        [roiMask, roiName] = load_roi_union(o.MaskDir, sub, ses, size(valid));
    end

    % ---- write ---------------------------------------------------------
    stem = fullfile(outDir, sprintf('%s_%s_run-%s', sub, ses, runLab));
    write_map(info, devImg,  [stem '_nbr-dev']);
    write_map(info, labImg,  [stem '_clust-label']);
    write_map(info, sizeImg, [stem '_clust-size']);
    write_map(info, oriImg,  [stem '_clust-ori']);

    M = sform_of(info);
    big = find(sz >= o.MinSize);
    for k = big(:)'
        mem = vi(lab == k);
        [cx, cy, cz] = ind2sub(size(valid), mem);
        ctr = [mean(cx) mean(cy) mean(cz)];
        w   = M * [ctr(:) - 1; 1];                 % sform is 0-based
        nRoi = 0;
        if ~isempty(roiMask), nRoi = sum(roiMask(mem)); end
        cRows = [cRows; {o.Trial, sub, ses, runLab, k, sz(k), ...
                         sz(k) * prod(double(info.PixelDimensions(1:3))), ...
                         oriC(k), Rin(k), median(devV(lab == k), 'omitnan'), ...
                         ctr(1), ctr(2), ctr(3), w(1), w(2), w(3), ...
                         nRoi, roiName}];                           %#ok<AGROW>
    end

    nRoiVox = 0; nRoiValid = 0; devRoi = NaN; maxRoi = 0;
    if ~isempty(roiMask)
        nRoiVox   = sum(roiMask(:));
        inRoi     = roiMask(vi);
        nRoiValid = sum(inRoi);
        if nRoiValid > 0
            devRoi = median(devV(inRoi), 'omitnan');
            maxRoi = max(sz(lab(inRoi)));
        end
    end

    fracBig = sz(1) / nv;                      % sz is sorted descending
    sRows = [sRows; {o.Trial, sub, ses, runLab, nRunsUsed, o.AngleTol, nv, ...
                     numel(ia), devObs, devNullMean, devP, ...
                     sz(1), p95, fracBig, Rin(1), sum(sz >= o.MinSize), ...
                     nRoiVox, nRoiValid, devRoi, maxRoi, roiName}];  %#ok<AGROW>

    fprintf('%-11s %-8s %7d %7.2f %8.2f %8.3f %9d %7.0f %6.3f %6.2f\n', ...
            sub, ses, nv, devObs, devNullMean, devP, sz(1), p95, fracBig, Rin(1));
    if fracBig > 0.25
        fprintf(2, ['  ^ the largest cluster holds %.0f%% of all valid ' ...
            'voxels. The field is PERCOLATING: that is one connected chain, ' ...
            'not a patch, and its extent is not a statistic. Check ' ...
            'max_clust_R (%.2f) - a real patch is near 1, a chain is near 0. ' ...
            'Smoothed functionals are the usual cause; lowering AngleTol is ' ...
            'the other lever. dev_deg stays readable either way.\n'], ...
            100 * fracBig, Rin(1));
    end
end

if ~isempty(cRows), append_clusters(csvPath, cRows); end
if ~isempty(sRows), append_summary(sumPath, sRows);  end
fprintf('\n%d clusters >= %d vox -> %s\n', size(cRows, 1), o.MinSize, csvPath);
fprintf('%d summary rows -> %s\n', size(sRows, 1), sumPath);
fprintf('maps -> %s\n', outDir);
fprintf('\nchance deviation is %g deg; smaller is smoother than chance.\n', 90/xf);

end

% ======================================================================
%  orientation field
% ======================================================================

function [u, v, valid, info, nRunsUsed] = orientation_field(gdir, event, xf, runSel)
% Unit orientation components. runSel 0 = across-run mean vector.
% Recovers (beta_cos, beta_sin) exactly from the two images GridCAT writes:
%   beta_cos = A*cos(x*ori*pi/180), beta_sin = A*sin(x*ori*pi/180)
% because ori = atan2(beta_sin,beta_cos)/x and A = hypot(beta_sin,beta_cos).
u = []; v = []; valid = []; info = []; nRunsUsed = 0;

f = dir(fullfile(gdir, sprintf('voxelwiseGridOri_%s_run*_deg.nii', event)));
bc = []; bs = []; nk = 0;
for k = 1:numel(f)
    n = regexp(f(k).name, '_run(\d+)_deg\.nii$', 'tokens', 'once');
    if isempty(n), continue; end
    r = str2double(n{1});
    if runSel ~= 0 && r ~= runSel, continue; end
    aPath = fullfile(gdir, sprintf('voxelwiseAmplitude_%s_run%d.nii', event, r));
    if ~isfile(aPath), continue; end
    if isempty(info), info = niftiinfo(fullfile(gdir, f(k).name)); end
    ori = double(niftiread(fullfile(gdir, f(k).name)));
    A   = double(niftiread(aPath));
    ang = xf * ori * pi / 180;
    nk  = nk + 1;
    bc(:,:,:,nk) = A .* cos(ang);                 %#ok<AGROW>
    bs(:,:,:,nk) = A .* sin(ang);                 %#ok<AGROW>
end
if nk == 0, return; end
nRunsUsed = nk;

% mean vector over the runs present, counting only voxels with an estimate.
% A == 0 is the atan2(0,0) trap: a zero-variance voxel reports ori 0 with no
% signal behind it, and must not be allowed to vote.
A  = hypot(bc, bs);
ok = isfinite(bc) & isfinite(bs) & A > 0;
bc(~ok) = 0; bs(~ok) = 0;
mc = sum(bc, 4); ms = sum(bs, 4);
cntr = sum(ok, 4);
mc(cntr < 1) = NaN; ms(cntr < 1) = NaN;

Am    = hypot(mc, ms);
valid = isfinite(Am) & Am > 0;
u = mc ./ Am; v = ms ./ Am;
u(~valid) = NaN; v(~valid) = NaN;
end

% ======================================================================
%  adjacency and components
% ======================================================================

function [a, b] = face_pairs(valid)
% Index pairs of face-adjacent valid voxels, in 1..nnz(valid) numbering.
% Face adjacency only: diagonals would make a 2-voxel-wide structure behave
% like a 3-voxel-wide one.
sz  = size(valid);
idx = zeros(sz);
vi  = find(valid);
idx(vi) = 1:numel(vi);
a = []; b = [];
for dim = 1:3
    if sz(dim) < 2, continue; end
    s1 = {':', ':', ':'};  s2 = s1;
    s1{dim} = 1:sz(dim)-1;
    s2{dim} = 2:sz(dim);
    A = idx(s1{:});  B = idx(s2{:});
    k = A > 0 & B > 0;
    a = [a; A(k)];  b = [b; B(k)];                %#ok<AGROW>
end
end

function [lab, sz] = components(a, b, nv)
% Connected components of "voxels that AGREE", which is a general graph and
% not a binary image - so not bwconncomp, and no Image Processing Toolbox.
% graph/conncomp are base MATLAB and are compiled, which matters: this runs
% once per permutation.
lab = conncomp(graph(a, b, [], nv))';
sz  = accumarray(lab, 1);
end

% ======================================================================
%  ROI, output, small helpers
% ======================================================================

function [U, name] = load_roi_union(maskDir, sub, ses, sz)
% Union of this session's ROI masks. The GLM's own mask.nii is SPM's
% implicit analysis mask, not an ROI, and is never read here.
U = []; name = '';
m = dir(fullfile(maskDir, sprintf('r%s_%s_hemi-*_mask.nii', sub, ses)));
if isempty(m), return; end
U = false(sz);
parts = {};
for k = 1:numel(m)
    X = niftiread(fullfile(m(k).folder, m(k).name));
    if numel(X) ~= prod(sz)
        error(['ROI mask and the maps are on different grids:\n' ...
               '  %s (%d voxels)\n  map (%d voxels)'], ...
               m(k).name, numel(X), prod(sz));
    end
    U = U | reshape(X > 0.5, sz);
    parts{end+1} = regexprep(m(k).name, ...
        '^r.*_hemi-([^_]+)_label-[^_]+_mask\.nii$', '$1');          %#ok<AGROW>
end
name = strjoin(unique(parts, 'stable'), '+');
end

function write_map(refInfo, X, stemNoExt)
info = refInfo;
info.Datatype              = 'single';
info.BitsPerPixel          = 32;
info.ImageSize             = size(X);
info.Description           = 'grid_clusters';
info.MultiplicativeScaling = 1;
info.AdditiveOffset        = 0;
niftiwrite(single(X), stemNoExt, info);
end

function M = sform_of(info)
r = info.raw;
M = double([r.srow_x; r.srow_y; r.srow_z; 0 0 0 1]);
if ~any(M(1:3, 1:3) ~= 0)
    error('This image has no sform; cannot place it in world space.');
end
end

function q = prctile_(x, pct)
% Linear-interpolated percentile, so the Statistics Toolbox is not needed.
x = sort(x(:));
n = numel(x);
if n == 0, q = NaN; return; end
if n == 1, q = x; return; end
h = (n - 1) * pct / 100 + 1;
lo = floor(h); hi = ceil(h);
q = x(lo) + (h - lo) * (x(hi) - x(lo));
end

function append_clusters(path, rows)
hdr = ['trial,subject,session,run,cluster,n_vox,vol_mm3,ori_deg,R_within,' ...
       'med_nbr_dev,cx,cy,cz,world_x,world_y,world_z,n_in_roi,roi'];
newFile = ~isfile(path);
fid = fopen(path, 'a');
if fid < 0, error('Cannot write CSV: %s', path); end
if newFile, fprintf(fid, '%s\n', hdr); end
for i = 1:size(rows, 1)
    fprintf(fid, '%s,%s,%s,%s,%d,%d,%.1f,%.3f,%.4f,%.3f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%d,%s\n', rows{i,:});
end
fclose(fid);
end

function append_summary(path, rows)
hdr = ['trial,subject,session,run,n_runs,angle_tol,n_valid,n_pairs,' ...
       'dev_deg,dev_null_deg,dev_p,max_clust,max_null_p95,frac_largest,' ...
       'max_clust_R,n_clust_min,n_roi,n_roi_valid,roi_dev_deg,' ...
       'roi_max_clust,roi'];
newFile = ~isfile(path);
fid = fopen(path, 'a');
if fid < 0, error('Cannot write CSV: %s', path); end
if newFile, fprintf(fid, '%s\n', hdr); end
for i = 1:size(rows, 1)
    fprintf(fid, '%s,%s,%s,%s,%d,%g,%d,%d,%.4f,%.4f,%.5f,%d,%.1f,%.4f,%.4f,%d,%d,%d,%.4f,%d,%s\n', rows{i,:});
end
fclose(fid);
end

function s = pick(c, a, b)
if c, s = a; else, s = b; end
end
