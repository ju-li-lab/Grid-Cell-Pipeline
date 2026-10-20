function info = nifti_info_via_fsl(niiFile)
% NIFTI_INFO_VIA_FSL  Get basic NIfTI info using fslinfo.
%   Returns struct with field 'nvols' (number of volumes).

[status, result] = system(sprintf('fslnvols %s', niiFile));
assert(status == 0, 'fslnvols failed for %s: %s', niiFile, result);
info.nvols = str2double(strtrim(result));
end
