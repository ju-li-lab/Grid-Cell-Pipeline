function coreg_reslice_rois(refEPI, t2w, roiFiles, CFG)
% Coregister: source = T2w, ref = mean EPI
% Reslice ROIs using the estimated transform (NN interp for masks).
% CFG is required for coregistration parameters
matlabbatch = {};

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {refEPI};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {t2w};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = roiFiles(:);
% Options from spm_options.m. roptions.interp defaults to 0 there
% (nearest neighbour), which is the correct choice for label masks.
roptions        = CFG.spm.coreg.roptions;
roptions.prefix = CFG.ROI_PREFIX;            % config.sh owns the prefix

matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions = CFG.spm.coreg.eoptions;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions = roptions;

spm_jobman('run', matlabbatch);
end
