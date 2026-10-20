function acqFile = resolve_acqparams(topupPrefix, funcDir, CFG)
% RESOLVE_ACQPARAMS  Locate the acqparams file applytopup needs.
%
% Prefers the one sitting next to the topup output, since that is what topup
% itself was run with. If topup output produced elsewhere came without one, a
% matching file is written into func/ from TOPUP_PE_DIR_* and TOPUP_READOUT_SEC.

acqFile = fullfile(fileparts(topupPrefix), 'topup_acqparams.txt');
if isfile(acqFile)
    return;
end

acqFile = fullfile(funcDir, 'topup_acqparams.txt');
if ~isfile(acqFile)
    fprintf('  No acqparams file next to the topup output — writing %s from config.\n', acqFile);
    write_acqparams(acqFile, CFG);
end
end
