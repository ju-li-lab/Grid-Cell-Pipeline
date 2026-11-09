function smooth_bold_files(funcDir, CFG)
% SMOOTH_BOLD_FILES  Apply Gaussian smoothing to preprocessed BOLD files.
%
% Finds the fully preprocessed BOLD files in funcDir, excluding mean images,
% and smooths them using SPM's smooth module.
% Output files get SMOOTH_PREFIX prepended (e.g. ua -> sua).
%
% WHICH FILES
% -----------
% CFG.BOLD_PREFIX, not RESLICE_PREFIX: the series to smooth is the one both
% realignment and slice timing have been through, and which letter comes
% first depends on SLICETIME_ORDER ('ua' or 'au'). The prefix is then
% followed by 'sub-', because a BIDS name starts there and matching on the
% prefix alone picks up anything that happens to begin with the same letter.
% SPM's FieldMap toolbox writes one such file per task - a single-volume
% unwarped sample named u<the EPI reference it used> - and an unanchored
% '^u.*_bold.*\.nii$' smoothed those too, reporting twice as many runs as
% the session has.

prefix = CFG.BOLD_PREFIX;  % 'u', 'ua' or 'au'

% Find preprocessed BOLD files but NOT mean images
d = dir(funcDir);
boldFiles = {};
for i = 1:numel(d)
    if d(i).isdir, continue; end
    % <prefix>sub-..._bold....nii, and never a mean image (mean*, <prefix>mean*)
    if ~isempty(regexp(d(i).name, ['^' regexptranslate('escape', prefix) 'sub-.*_bold.*\.nii$'], 'once')) && ...
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
