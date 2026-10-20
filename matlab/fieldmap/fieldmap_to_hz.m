function hzFile = fieldmap_to_hz(fieldmapFile, funcDir, CFG)
% FIELDMAP_TO_HZ  Write a copy of the fieldmap in Hz for SPM's FieldMap toolbox.
%
% fsl_prepare_fieldmap writes rad/s; SPM's "Precalculated FieldMap" input is
% documented as Hz, so rad/s is divided by 2*pi. A fieldmap already in Hz is
% copied through unchanged, so the rest of the code has one path to follow.
%
% Doing the unit conversion here — and letting SPM do the Hz -> voxel-shift
% conversion with its own internal scaling (from TOTAL_READOUT_MS and
% BLIP_DIRECTION) — keeps this code out of the business of guessing what
% scaling SPM's VDM files use.

V = spm_vol(fieldmapFile);
if numel(V) > 1, V = V(1); end
dat = spm_read_vols(V);

switch CFG.FIELDMAP_UNITS
    case {'rad/s','rads','rad_per_s'}
        dat = dat / (2*pi);
        fprintf('  Converted fieldmap rad/s -> Hz (divided by 2*pi).\n');
    case 'hz'
        fprintf('  Fieldmap is already in Hz; no conversion.\n');
    otherwise
        error('Unknown FIELDMAP_UNITS: %s (must be rad/s or Hz)', CFG.FIELDMAP_UNITS);
end

fprintf('  Fieldmap range: %.1f to %.1f Hz\n', min(dat(:)), max(dat(:)));

hzFile = fullfile(funcDir, 'fieldmap_hz.nii');
Vout = V;
Vout.fname = hzFile;
Vout.dt = [spm_type('float32') 0];
Vout.pinfo = [1;0;0];
Vout.descrip = 'B0 fieldmap in Hz (for SPM FieldMap)';
spm_write_vol(Vout, dat);
end
