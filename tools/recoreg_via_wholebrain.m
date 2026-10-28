function recoreg_via_wholebrain(bidsRoot, derivRoot, sub, ses, roiLabel)
% RECOREG_VIA_WHOLEBRAIN  Redo one session's T2w->EPI coregistration, going
% through the whole-brain resting-state EPI instead of straight to the task
% EPI.
%
%   recoreg_via_wholebrain(bidsRoot, derivRoot, 'sub-05s17', 'ses-01', 'princeErC')
%
% WHY
%   The T2w slab (175x175x80 mm) and the task EPI slab (192x192x104 mm) only
%   share a small block of tissue, and normalised mutual information over a
%   small overlap has a degenerate optimum: driving the images apart shrinks
%   the overlap to a near-uniform sliver that can score better than the true
%   alignment. That is how sub-05 ses-01 ended up translated ~141 mm, with
%   both ROI masks resliced to zero voxels.
%
%   The resting-state EPI covers 216x216x158 mm, so it overlaps both images
%   generously. Stangl et al. (2018) used a whole-brain EPI for exactly this
%   reason - their own methods call the direct route "non-trivial".
%
% THE MASKS ARE RESAMPLED ONCE, NOT TWICE
%   SPM separates estimate from reslice. Both steps below are ESTIMATE only,
%   which just updates the orientation matrix in each header. The single
%   write at the end applies the composed transform. Two rigid transforms
%   compose into one rigid transform, so this costs no extra interpolation
%   over the one-step route.
%
% NON-DESTRUCTIVE
%   Everything happens in <deriv>/<sub>/<ses>/anat/recoreg/ on fresh copies
%   taken from the raw BIDS tree. Nothing the pipeline wrote is touched. The
%   function prints the one command that installs the result if you like it.
%
% This is a scratch tool for trying the idea on one session. If it works it
% belongs in matlab/steps/, behind a config flag.

if nargin < 5, roiLabel = 'princeErC'; end

rawAnat = fullfile(bidsRoot,  sub, ses, 'anat');
rawFunc = fullfile(bidsRoot,  sub, ses, 'func');
derAnat = fullfile(derivRoot, sub, ses, 'anat');
derFunc = fullfile(derivRoot, sub, ses, 'func');
work    = fullfile(derAnat, 'recoreg');
if ~isfolder(work), mkdir(work); end

% ---- inputs -------------------------------------------------------------
t2wSrc  = pick_one(rawAnat, sprintf('%s_%s*_T2w.nii*', sub, ses), 'T2w');
maskSrc = dir(fullfile(rawAnat, sprintf('%s_%s*label-%s_mask.nii*', sub, ses, roiLabel)));
assert(~isempty(maskSrc), 'No %s masks in %s', roiLabel, rawAnat);
restSrc = pick_one(rawFunc, sprintf('%s_%s_task-rest*_bold.nii*', sub, ses), ...
                   'whole-brain (rest) EPI');
meanSrc = pick_one(derFunc, 'meanu*.nii', 'session mean EPI');

fprintf('\n=== recoreg %s %s (%s) ===\n', sub, ses, roiLabel);
fprintf('  T2w          : %s\n', shortname(t2wSrc));
fprintf('  whole-brain  : %s\n', shortname(restSrc));
fprintf('  session mean : %s\n', shortname(meanSrc));

% ---- working copies (SPM rewrites headers, so never the originals) ------
t2w   = stage(t2wSrc, work);
masks = arrayfun(@(m) stage(fullfile(m.folder, m.name), work), maskSrc, ...
                 'UniformOutput', false);
wb    = first_volume(stage(restSrc, work), fullfile(work, 'wholebrain_epi.nii'));

fprintf('\n  before:  T2w is %.0f mm from the session mean\n', ...
        com_dist(t2w, meanSrc));
for k = 1:numel(masks)
    old = fullfile(derAnat, ['r' spm_file(masks{k}, 'filename')]);
    if isfile(old)
        fprintf('           %-38s %5d voxels (pipeline)\n', ...
                spm_file(old, 'filename'), nvox(old));
    end
end

% ---- SPM's own defaults; the failure is the geometry, not these ---------
eo = struct('cost_fun', 'nmi', 'sep', [4 2 1], ...
            'tol', [0.02 0.02 0.02 0.001 0.001 0.001 ...
                    0.01 0.01 0.01 0.001 0.001 0.001], ...
            'fwhm', [7 7]);

% ---- step 1: T2w (+ masks) -> whole-brain EPI    ESTIMATE ONLY ----------
run_batch(estimate_batch(wb, t2w, masks, eo));

% ---- step 2: whole-brain (+ T2w + masks) -> session mean   ESTIMATE -----
run_batch(estimate_batch(meanSrc, wb, [{t2w}; masks(:)], eo));

% ---- step 3: the single reslice ----------------------------------------
mb = {};
mb{1}.spm.spatial.coreg.write.ref      = {meanSrc};
mb{1}.spm.spatial.coreg.write.source   = masks(:);
mb{1}.spm.spatial.coreg.write.roptions = struct( ...
    'interp', 0, 'wrap', [0 0 0], 'mask', 0, 'prefix', 'r');   % NN for labels
run_batch(mb);

% a trilinear T2w as well, purely to look at
mb = {};
mb{1}.spm.spatial.coreg.write.ref      = {meanSrc};
mb{1}.spm.spatial.coreg.write.source   = {t2w};
mb{1}.spm.spatial.coreg.write.roptions = struct( ...
    'interp', 4, 'wrap', [0 0 0], 'mask', 0, 'prefix', 'r');
run_batch(mb);

% ---- report -------------------------------------------------------------
fprintf('\n  after:   T2w is %.0f mm from the session mean\n', com_dist(t2w, meanSrc));
ok = true;
for k = 1:numel(masks)
    out = fullfile(work, ['r' spm_file(masks{k}, 'filename')]);
    n   = nvox(out);
    fprintf('           %-38s %5d voxels\n', spm_file(out, 'filename'), n);
    if n < 5, ok = false; end
end

if ok
    fprintf(['\n  Looks sane. To use it, overwrite the pipeline''s masks:\n' ...
             '    cp %s/r*_mask.nii %s/\n'], work, derAnat);
else
    fprintf(['\n  Still near-empty. Open %s/r%s in a viewer on top of the\n' ...
             '  session mean before trusting any of this.\n'], ...
             work, spm_file(t2w, 'filename'));
end
end


% =========================================================================
function mb = estimate_batch(ref, src, other, eo)
mb = {};
mb{1}.spm.spatial.coreg.estimate.ref      = {ref};
mb{1}.spm.spatial.coreg.estimate.source   = {src};
mb{1}.spm.spatial.coreg.estimate.other    = other(:);
mb{1}.spm.spatial.coreg.estimate.eoptions = eo;
end

function run_batch(mb)
spm_jobman('run', mb);
end

function f = pick_one(dirPath, pat, what)
d = dir(fullfile(dirPath, pat));
d = d(~[d.isdir]);
assert(~isempty(d), 'No %s matching %s in %s', what, pat, dirPath);
f = fullfile(d(1).folder, d(1).name);
if numel(d) > 1
    fprintf('  (%d files match %s, using %s)\n', numel(d), pat, d(1).name);
end
end

function out = stage(src, destDir)
% Copy into the working directory, decompressing if needed. SPM cannot read
% .gz and will rewrite whatever header it is given.
[~, base, ext] = fileparts(src);
if strcmpi(ext, '.gz')
    produced = gunzip(src, destDir);
    out = produced{1};
else
    out = fullfile(destDir, [base ext]);
    copyfile(src, out);
end
fileattrib(out, '+w');
end

function out = first_volume(src4d, out)
V = spm_vol([src4d ',1']);
Y = spm_read_vols(V);
% .private still points at the 4D source; leaving it makes spm_write_vol
% write back into that file instead of the new one.
if isfield(V, 'private'), V = rmfield(V, 'private'); end
V.fname = out; V.n = [1 1];
spm_write_vol(V, Y);
end

function n = nvox(f)
if ~isfile(f), n = -1; return; end
Y = spm_read_vols(spm_vol(f));
n = sum(Y(:) > 0.5);
end

function d = com_dist(a, b)
d = norm(com(a) - com(b));
end

function c = com(f)
V = spm_vol(f); if numel(V) > 1, V = V(1); end
Y = spm_read_vols(V);
Y(~isfinite(Y) | Y < 0) = 0;
[i, j, k] = ndgrid(1:size(Y,1), 1:size(Y,2), 1:size(Y,3));
w = Y(:); s = sum(w);
if s == 0, c = [0 0 0]; return; end
vox = [sum(i(:).*w) sum(j(:).*w) sum(k(:).*w)] / s;
c = (V.mat * [vox 1]')';
c = c(1:3);
end

function s = shortname(p)
[d, n, e] = fileparts(p);
[~, dn] = fileparts(d);
s = fullfile(dn, [n e]);
end
