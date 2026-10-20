function vdm_out = convert_topup_to_vdm(topupPrefix, epiRef, funcDir, CFG, taskLabel)
% CONVERT_TOPUP_TO_VDM  Convert FSL topup field to SPM-compatible VDM.
%
% The topup field (in Hz) is converted to a voxel displacement map (in mm)
% using the formula: VDM = fieldmap_Hz * totalReadoutTime * voxelSize_PE
%
% The VDM is then written as a NIfTI in funcDir for use with SPM Realign&Unwarp.

fieldFile = sprintf('%s_field.nii.gz', topupPrefix);
if ~isfile(fieldFile)
    % Topup output made outside this pipeline may not be compressed
    fieldFile = sprintf('%s_field.nii', topupPrefix);
end
assert(isfile(fieldFile), 'Topup field not found: %s_field.nii[.gz]', topupPrefix);

% Convert to .nii for SPM (FSL defaults to .nii.gz)
fieldNii = fullfile(funcDir, sprintf('topup_field_task-%s.nii', taskLabel));
fieldNiiGz = [fieldNii '.gz'];
cmd = sprintf('fslmaths %s %s -odt float', fieldFile, fieldNiiGz);
run_shell(cmd);
% Decompress .nii.gz -> .nii
if isfile(fieldNiiGz)
    gunzip(fieldNiiGz);
    delete(fieldNiiGz);
end
assert(isfile(fieldNii), 'fslmaths output not found: %s', fieldNii);

% Read the field with SPM
Vfield = spm_vol(fieldNii);
fieldData = spm_read_vols(Vfield);

% Read EPI reference to get geometry (first volume only)
% epiRef may contain ',1' frame specifier; if not, spm_vol on 4D returns
% an array of structs — take only the first element.
epiPath = regexp(epiRef, '^[^,]+', 'match', 'once');
Vepi = spm_vol(epiPath);
if numel(Vepi) > 1
    Vepi = Vepi(1);
end

% Determine PE direction index and voxel size along that axis
peDir = CFG.TOPUP_PE_DIR_BOLD;
switch peDir
    case {'x', 'x-'}
        peDim = 1;
    case {'y', 'y-'}
        peDim = 2;
    case {'z', 'z-'}
        peDim = 3;
    otherwise
        error('Unknown PE direction: %s', peDir);
end
voxSize_PE = abs(Vepi.mat(peDim, peDim));

% Convert the topup field (Hz) to a voxel displacement map.
%
% SPM's VDM is in VOXELS: spm_uw_apply shifts each voxel by the VDM value
% along the phase-encode direction, in voxel units. field_Hz * readout_s is
% already that number, so no voxel-size factor belongs here. An earlier
% version multiplied by voxSize_PE, which over-corrected the distortion by
% the PE voxel size (2-3x for a typical EPI).
readout_sec = CFG.TOPUP_READOUT_SEC;
vdmData = fieldData * readout_sec;

% If PE direction is negative, flip the sign.
%
% UNVERIFIED: topup's --fout sign is already defined relative to the first
% row of acqparams.txt, so this flip may double-count. If the unwarped
% images come out MORE distorted than the originals rather than less, try
% removing this branch - that is what a sign error looks like.
if endsWith(peDir, '-')
    vdmData = -vdmData;
end

% Report the peak shift, so an implausible VDM is visible in the job log.
% A sane EPI distortion is a few voxels; tens of voxels means the readout
% time or the field units are wrong.
fprintf('  VDM task-%s: peak shift %.2f voxels (= %.2f mm at %.2f mm PE voxel)\n', ...
    taskLabel, max(abs(vdmData(:))), max(abs(vdmData(:))) * voxSize_PE, voxSize_PE);

% Write VDM using field header (same geometry)
vdm_out = fullfile(funcDir, sprintf('vdm_task-%s.nii', taskLabel));
Vout = Vfield;
Vout.fname = vdm_out;
Vout.dt = [spm_type('float32') 0];
Vout.descrip = 'VDM from FSL topup';
spm_write_vol(Vout, vdmData);

% Clean up intermediate field file
if isfile(fieldNii), delete(fieldNii); end

end
