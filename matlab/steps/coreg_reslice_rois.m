function coreg_reslice_rois(refEPI, t2w, roiFiles, CFG, bridgeEPI)
% Coregister the T2w to the session mean and reslice it and the ROI masks
% into EPI space.
%
% Two routes, chosen by COREG_VIA_TASK in config.sh:
%
%   direct      T2w -> session mean                      (COREG_VIA_TASK empty)
%   bridged     T2w -> whole-brain EPI -> session mean   (COREG_VIA_TASK=rest)
%
% The bridged route exists because a T2w slab and a partial-volume task EPI
% overlap in only a small block of tissue. Normalised mutual information over
% a small overlap has a degenerate optimum: driving the images apart shrinks
% the overlap to a near-uniform sliver, which can score better than the true
% alignment. That is how sub-05 ses-01 ended up translated 141 mm with both
% masks resliced to zero voxels. A whole-brain EPI overlaps both images
% generously, so neither step is a hard registration.
%
% THE MASKS ARE RESAMPLED ONCE, NOT TWICE.
% Both coregistration calls below are ESTIMATE only - they just update the
% orientation matrix in each header and move no voxels. The single write at
% the end applies whatever transform has accumulated there. Two rigid
% transforms compose into one rigid transform, so the bridged route costs no
% interpolation over the direct one.
%
% Everything that must follow the T2w rides along as `other`, so the masks
% never get a transform of their own.

if nargin < 5, bridgeEPI = ''; end

rois   = roiFiles(:);
moving = [{t2w}; rois];        % reslice the T2w too, for looking at

if isempty(bridgeEPI)
    estimate(refEPI, t2w, rois, CFG);
else
    bridge = first_volume(bridgeEPI);
    fprintf('  bridging through %s\n', spm_file(bridge, 'filename'));
    estimate(bridge,  t2w,    rois,   CFG);   % T2w (+ masks) -> bridge
    estimate(refEPI,  bridge, moving, CFG);   % bridge (+ T2w + masks) -> mean
end

reslice(refEPI, moving, CFG);
check_rois(rois, CFG);
end


% -------------------------------------------------------------------------
function estimate(ref, source, other, CFG)
% Header-only. Nothing is resampled here.
if isempty(other), other = {''}; end
matlabbatch = {};
matlabbatch{1}.spm.spatial.coreg.estimate.ref      = {ref};
matlabbatch{1}.spm.spatial.coreg.estimate.source   = {source};
matlabbatch{1}.spm.spatial.coreg.estimate.other    = other(:);
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions = CFG.spm.coreg.eoptions;
spm_jobman('run', matlabbatch);
end


function reslice(ref, source, CFG)
% The one resampling. roptions.interp is 0 in spm_options.m (nearest
% neighbour), which is the correct choice for label masks.
roptions        = CFG.spm.coreg.roptions;
roptions.prefix = CFG.ROI_PREFIX;            % config.sh owns the prefix

matlabbatch = {};
matlabbatch{1}.spm.spatial.coreg.write.ref      = {ref};
matlabbatch{1}.spm.spatial.coreg.write.source   = source(:);
matlabbatch{1}.spm.spatial.coreg.write.roptions = roptions;
spm_jobman('run', matlabbatch);
end


function out = first_volume(nii4d)
% SPM wants one volume to register against, and the second estimate rewrites
% its header - which must not mean rewriting a frame of the raw 4D run.
V = spm_vol(nii4d);
if numel(V) == 1, out = nii4d; return; end

out = spm_file(nii4d, 'prefix', 'vol1_');
matlabbatch = {};
matlabbatch{1}.spm.util.imcalc.input          = {sprintf('%s,1', nii4d)};
matlabbatch{1}.spm.util.imcalc.output         = spm_file(out, 'filename');
matlabbatch{1}.spm.util.imcalc.outdir         = {spm_file(out, 'path')};
matlabbatch{1}.spm.util.imcalc.expression     = 'i1';
matlabbatch{1}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype  = 16;
spm_jobman('run', matlabbatch);
end


function check_rois(rois, CFG)
% A coregistration that goes wrong does not raise an error. It writes a mask
% that misses the brain, GridCAT then finds nothing in it, and the session
% quietly vanishes from the analysis with no record of why. Stopping here
% puts the failure where the cause is.
minVox = CFG.COREG_MIN_ROI_VOXELS;
for k = 1:numel(rois)
    out = spm_file(rois{k}, 'prefix', CFG.ROI_PREFIX);
    n   = 0;
    if isfile(out)
        Y = spm_read_vols(spm_vol(out));
        n = sum(Y(:) > 0.5);
    end
    fprintf('  %-46s %5d voxels\n', spm_file(out, 'filename'), n);

    if minVox > 0 && n < minVox
        error('coreg:EmptyROI', [ ...
            '%s came back with %d voxels (COREG_MIN_ROI_VOXELS is %d).\n\n' ...
            'The coregistration has put the T2w somewhere the EPI does not ' ...
            'cover, so the mask missed the brain entirely.\n' ...
            'Look at\n  %s\non top of the session mean.\n\n' ...
            'If the T2w and the task EPI are both partial slabs, set ' ...
            'COREG_VIA_TASK in config.sh to a whole-brain run (usually ' ...
            '"rest") and run --stages coreg again.'], ...
            spm_file(out, 'filename'), n, minVox, out);
    end
end
end
