function smooth_bold_files(funcDir, CFG)
% SMOOTH_BOLD_FILES  Apply Gaussian smoothing to preprocessed BOLD files.
%
% Finds all resliced BOLD files (u*_bold*.nii) in funcDir, excluding mean
% images, and smooths them using SPM's smooth module.
% Output files get SMOOTH_PREFIX prepended (e.g. u -> su).

prefix = CFG.RESLICE_PREFIX;  % typically 'u'

% Find preprocessed BOLD files but NOT mean images
d = dir(funcDir);
boldFiles = {};
for i = 1:numel(d)
    if d(i).isdir, continue; end
    % Match files starting with the reslice prefix, containing _bold, ending .nii
    % but exclude mean images (meanu*, mean*)
    if ~isempty(regexp(d(i).name, ['^' prefix '.*_bold.*\.nii$'], 'once')) && ...
       isempty(regexp(d(i).name, '^mean', 'once'))
        boldFiles{end+1} = fullfile(funcDir, d(i).name); %#ok<AGROW>
    end
end

if isempty(boldFiles)
    warning('smooth_bold_files:NoFiles', 'No preprocessed BOLD files found to smooth in %s', funcDir);
    return;
end

fprintf('  Found %d BOLD files to smooth.\n', numel(boldFiles));

% Ensure isotropic kernel
fwhm = CFG.SMOOTH_FWHM;
if isscalar(fwhm)
    fwhm = [fwhm fwhm fwhm];
end

matlabbatch = {};
matlabbatch{1}.spm.spatial.smooth.data   = boldFiles(:);
matlabbatch{1}.spm.spatial.smooth.fwhm   = fwhm;
matlabbatch{1}.spm.spatial.smooth.dtype  = CFG.spm.smooth.dtype;
matlabbatch{1}.spm.spatial.smooth.im     = CFG.spm.smooth.im;
matlabbatch{1}.spm.spatial.smooth.prefix = CFG.SMOOTH_PREFIX;

spm_jobman('run', matlabbatch);
fprintf('  Smoothing complete. Output prefix: %s\n', CFG.SMOOTH_PREFIX);
end

% ================= PRECALCULATED FIELDMAP HELPERS =================
