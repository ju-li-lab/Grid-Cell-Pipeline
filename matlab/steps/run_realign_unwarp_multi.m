function run_realign_unwarp_multi(taskVols, taskVDM, CFG)
% ONE Realign&Unwarp batch, multiple data blocks.
% taskVols: cell array, each element is cellstr of scans for that run
% taskVDM : cell array, each element is full path to vdm for that run
% CFG     : configuration struct with realign/unwarp parameters

n = numel(taskVols);
assert(numel(taskVDM) == n, 'taskVols and taskVDM must have same length.');

matlabbatch = {};
for i = 1:n
    matlabbatch{1}.spm.spatial.realignunwarp.data(i).scans  = taskVols{i};
    matlabbatch{1}.spm.spatial.realignunwarp.data(i).pmscan = {taskVDM{i}};
end

% Options come straight from spm_options.m. Add any field SPM accepts
% there and it arrives here without a change to this file.
uwroptions        = CFG.spm.realignunwarp.uwroptions;
uwroptions.prefix = CFG.RESLICE_PREFIX;      % config.sh owns the prefix

matlabbatch{1}.spm.spatial.realignunwarp.eoptions   = CFG.spm.realignunwarp.eoptions;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions = CFG.spm.realignunwarp.uweoptions;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions = uwroptions;

spm_jobman('run', matlabbatch);
end
