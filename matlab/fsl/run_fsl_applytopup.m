function outputFile = run_fsl_applytopup(boldFile, topupPrefix, acqparamsFile, outputFile, imainIndex, CFG)
% RUN_FSL_APPLYTOPUP  Apply topup distortion correction to a 4D BOLD file.
%
% Inputs:
%   boldFile      - Input 4D BOLD NIfTI
%   topupPrefix   - Prefix from run_fsl_topup (contains fieldcoefs, movpar)
%   acqparamsFile - acqparams.txt topup was run with (see resolve_acqparams)
%   outputFile    - Desired output path (must end in .nii)
%   imainIndex    - Index into acqparams.txt for this image (1 = forward PE)
%   CFG           - Config struct with TOPUP_INTERP
%
% Returns:
%   outputFile  - Actual path to the corrected file (.nii, decompressed)
%
% Notes:
%   - FSL --out expects a basename (no .nii extension). We strip it.
%   - FSL always writes .nii.gz output. We decompress to .nii for SPM.
%   - We use --method=jac (Jacobian modulation) for single-image correction.
%     The alternative --method=lsr (least-squares) requires BOTH forward and
%     reverse PE images passed together, which doesn't apply here.

assert(isfile(acqparamsFile), 'acqparams file not found: %s', acqparamsFile);

% FSL --out expects a basename without extension
outBase = regexprep(outputFile, '\.nii(\.gz)?$', '');

cmd = sprintf('applytopup --imain=%s --topup=%s --datain=%s --inindex=%d --method=jac --interp=%s --out=%s', ...
    boldFile, topupPrefix, acqparamsFile, imainIndex, CFG.TOPUP_INTERP, outBase);

fprintf('  Applying topup to %s:\n    %s\n', boldFile, cmd);
run_shell(cmd);

% FSL writes <basename>.nii.gz — decompress to .nii for SPM
outGz = [outBase '.nii.gz'];
outputFile = [outBase '.nii'];

if isfile(outGz)
    fprintf('  Decompressing applytopup output for SPM...\n');
    gunzip(outGz);
    delete(outGz);
end

assert(isfile(outputFile), 'applytopup output not found: %s (checked %s too)', outputFile, outGz);
fprintf('  applytopup complete -> %s\n', outputFile);
end
