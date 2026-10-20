function write_acqparams(acqFile, CFG)
% WRITE_ACQPARAMS  Write the two-line FSL acqparams file from the config.
% Line 1 = forward (BOLD) phase-encode direction, line 2 = reverse.

fid = fopen(acqFile, 'w');
assert(fid > 0, 'Cannot create acqparams file: %s', acqFile);

fwdVec = pe_dir_to_vector(CFG.TOPUP_PE_DIR_BOLD);
revVec = pe_dir_to_vector(CFG.TOPUP_PE_DIR_REVERSE);
readout = CFG.TOPUP_READOUT_SEC;

fprintf(fid, '%d %d %d %.6f\n', fwdVec(1), fwdVec(2), fwdVec(3), readout);
fprintf(fid, '%d %d %d %.6f\n', revVec(1), revVec(2), revVec(3), readout);
fclose(fid);
end
