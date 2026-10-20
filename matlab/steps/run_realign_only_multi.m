function run_realign_only_multi(taskVols, CFG)
% Realign only (no unwarp/VDM).
% taskVols: cell array, each element is cellstr of scans for that run
% CFG     : configuration struct with realign parameters

n = numel(taskVols);

matlabbatch = {};

% Data specification
for i = 1:n
    matlabbatch{1}.spm.spatial.realign.estwrite.data{i} = taskVols{i};
end

% Options from spm_options.m. Note SPM names these differently from
% Realign & Unwarp: 'interp'/'wrap' here, 'einterp'/'ewrap' there.
roptions        = CFG.spm.realign.roptions;
roptions.prefix = CFG.RESLICE_PREFIX;        % config.sh owns the prefix

matlabbatch{1}.spm.spatial.realign.estwrite.eoptions = CFG.spm.realign.eoptions;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions = roptions;

spm_jobman('run', matlabbatch);
end
