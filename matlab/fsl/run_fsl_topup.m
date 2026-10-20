function run_fsl_topup(boldFile, reversePE_epi, topupPrefix, CFG)
% RUN_FSL_TOPUP  Estimate distortion field using FSL topup.
%
% Creates an acqparams.txt file, merges the first BOLD volume with the
% reverse-PE EPI, and runs FSL topup to estimate the fieldmap.
%
% Inputs:
%   boldFile      - Path to the first task's 4D BOLD (uses volume 1)
%   reversePE_epi - Path to the reverse phase-encode EPI in fmap/
%   topupPrefix   - Output prefix for topup results (e.g., funcDir/topup_results)
%   CFG           - Config struct with TOPUP_* fields

funcDir = fileparts(topupPrefix);

% --- Step 1: Extract first volume from BOLD as the "forward" EPI ---
forwardEPI = fullfile(funcDir, 'topup_forward_b0.nii.gz');
cmd = sprintf('fslroi %s %s 0 1', boldFile, forwardEPI);
fprintf('  Extracting first BOLD volume: %s\n', cmd);
run_shell(cmd);

% Handle reverse-PE EPI: might be 4D, take first volume
reverseInfo = nifti_info_via_fsl(reversePE_epi);
if reverseInfo.nvols > 1
    reverseB0 = fullfile(funcDir, 'topup_reverse_b0.nii.gz');
    cmd = sprintf('fslroi %s %s 0 1', reversePE_epi, reverseB0);
    fprintf('  Extracting first reverse-PE volume: %s\n', cmd);
    run_shell(cmd);
else
    reverseB0 = reversePE_epi;
end

% --- Step 2: Merge forward + reverse into a single 4D file ---
mergedFile = fullfile(funcDir, 'topup_merged_b0.nii.gz');
cmd = sprintf('fslmerge -t %s %s %s', mergedFile, forwardEPI, reverseB0);
fprintf('  Merging forward + reverse EPIs: %s\n', cmd);
run_shell(cmd);

% --- Step 3: Create acqparams.txt ---
acqparamsFile = fullfile(funcDir, 'topup_acqparams.txt');
write_acqparams(acqparamsFile, CFG);
fprintf('  Acqparams file: %s\n', acqparamsFile);

% --- Step 4: Run FSL topup ---
topupCmd = sprintf('topup --imain=%s --datain=%s --out=%s --fout=%s_field --iout=%s_corrected', ...
    mergedFile, acqparamsFile, topupPrefix, topupPrefix, topupPrefix);

if ~isempty(CFG.TOPUP_CONFIG)
    topupCmd = sprintf('%s --config=%s', topupCmd, CFG.TOPUP_CONFIG);
end

% Add verbose flag for logging
topupCmd = sprintf('%s --verbose', topupCmd);

fprintf('  Running FSL topup:\n    %s\n', topupCmd);
run_shell(topupCmd);
fprintf('  Topup estimation complete.\n');

% Clean up temporary files
if isfile(forwardEPI), delete(forwardEPI); end
if ~strcmp(reverseB0, reversePE_epi) && isfile(reverseB0), delete(reverseB0); end
if isfile(mergedFile), delete(mergedFile); end

end
