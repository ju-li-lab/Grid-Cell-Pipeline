function vdm_out = make_vdm(data, epiRef, funcDir, fmapDir, CFG, taskLabel)
% MAKE_VDM  Build one voxel displacement map with SPM's FieldMap toolbox.
%
%   vdm_out = make_vdm(data, epiRef, funcDir, fmapDir, CFG, taskLabel)
%
% DATA says which kind of input SPM is being given, and is one of:
%
%   struct('kind','presubphasemag', 'phase',<phasediff>, 'magnitude',<mag1>)
%       A GRE phase-difference image and its magnitude. SPM unwraps the
%       phase itself. This is PREPROC_MODE=realign_unwarp.
%
%   struct('kind','precalcfieldmap', 'fieldmap',<fieldmap_Hz>, 'magnitude',<mag>)
%       A ready-made fieldmap, already in Hz. SPM only converts it to a
%       VDM. This is PREPROC_MODE=precalc_fieldmap.
%
% Everything else is identical between the two, which is why they are one
% function: the GRE and precalculated routes were previously two
% near-verbatim copies, so any fieldmap fix had to be made twice and the
% hardcoded mflags existed in two places.
%
% The result is moved to funcDir under a deterministic name,
% vdm_task-<label>.nii, because vdmflags.prefix cannot be relied on -
% SPM's FieldMap toolbox sometimes ignores or overrides it. So instead of
% trusting the name, this snapshots the vdm* files before the run and diffs
% afterwards to find what SPM actually wrote.

% ---- Snapshot existing VDM candidates (before) ----
vdmSearchDirs = {funcDir, fmapDir, pwd};
before = list_vdm_candidates(vdmSearchDirs);

% ---- Build the FieldMap batch ----
matlabbatch = {};

switch data.kind
    case 'presubphasemag'
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase     = {data.phase};
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = {data.magnitude};
    case 'precalcfieldmap'
        % Echo times play no part here - there is no phase to unwrap - but
        % the batch still requires a complete defaults structure below.
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.precalcfieldmap.precalcfieldmap = {data.fieldmap};
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.precalcfieldmap.magfieldmap     = {data.magnitude};
    otherwise
        error('make_vdm: unknown data kind ''%s'' (expected presubphasemag or precalcfieldmap).', data.kind);
end

matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session.epi = {epiRef};

% Everything SPM's FieldMap toolbox takes, straight from spm_options.m.
% Add any field SPM accepts there and it arrives here unchanged.
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval = fieldmap_defaults(CFG);
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.vdmflags      = CFG.spm.fieldmap.vdmflags;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm      = CFG.spm.fieldmap.matchvdm;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = CFG.spm.fieldmap.writeunwarped;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat     = CFG.spm.fieldmap.matchanat;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname      = sprintf('task-%s', taskLabel);
% Avoids "incomplete module inputs" errors on some installs.
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat          = {''};

% ---- Run from funcDir, so outputs are more likely to land there ----
oldpwd  = pwd;
cleanup = onCleanup(@() cd(oldpwd)); %#ok<NASGU>
cd(funcDir);

spm_jobman('run', matlabbatch);

% ---- Snapshot after, diff to find the new VDM ----
after    = list_vdm_candidates([vdmSearchDirs, {pwd}]);
newFiles = diff_vdm_candidates(before, after);
assert(~isempty(newFiles), ...
    'VDM not found after the FieldMap run for task-%s. Check where FieldMap writes its output.', taskLabel);

newest = pick_newest_file(newFiles);   % if several appeared, take the newest

% ---- Move to a deterministic name in funcDir ----
vdm_out = fullfile(funcDir, sprintf('vdm_task-%s.nii', taskLabel));

[~, ~, ext] = fileparts(newest);
if strcmpi(ext, '.img')
    % Analyze .img/.hdr pair: move both, keep the .img name. SPM reads it.
    target_img = fullfile(funcDir, sprintf('vdm_task-%s.img', taskLabel));
    target_hdr = fullfile(funcDir, sprintf('vdm_task-%s.hdr', taskLabel));
    if isfile(target_img), delete(target_img); end
    if isfile(target_hdr), delete(target_hdr); end
    movefile(newest, target_img);
    hdr = strrep(newest, '.img', '.hdr');
    if isfile(hdr), movefile(hdr, target_hdr); end
    vdm_out = target_img;
else
    if isfile(vdm_out), delete(vdm_out); end
    movefile(newest, vdm_out);
end

fprintf('VDM (task-%s) -> %s\n', taskLabel, vdm_out);

end


function d = fieldmap_defaults(CFG)
% SPM FieldMap defaults, with the brain-mask template filled in.
%
% spm_options.m leaves mflags.template empty so that reading the settings
% does not require SPM to be on the path. SPM's own T1 is the default.
d = CFG.spm.fieldmap.defaultsval;
if isempty(d.mflags.template) || isempty(d.mflags.template{1})
    d.mflags.template = {fullfile(spm('Dir'), 'toolbox', 'FieldMap', 'T1.nii')};
end
end


function files = list_vdm_candidates(folders)
% LIST_VDM_CANDIDATES  Snapshot the vdm* files across one or more folders.
%
% Returns a struct array with .path, .datenum and .bytes. Recording the
% timestamp and size (not just the path) means a VDM that SPM overwrote in
% place still registers as new output — which happens when a stale vdm5_* from
% an earlier manual run is already sitting in fmap/.

if ischar(folders), folders = {folders}; end
files = struct('path', {}, 'datenum', {}, 'bytes', {});
seen = {};

for f = 1:numel(folders)
    folder = folders{f};
    if isempty(folder) || ~isfolder(folder), continue; end
    d = dir(folder);
    for i = 1:numel(d)
        if d(i).isdir, continue; end
        n = d(i).name;

        % FieldMap typically uses vdm*.nii or vdm*.img/.hdr
        if isempty(regexp(n,'^vdm.*\.(nii|img)$','once')), continue; end

        p = fullfile(folder, n);
        if ismember(p, seen), continue; end
        seen{end+1} = p; %#ok<AGROW>
        files(end+1) = struct('path', p, 'datenum', d(i).datenum, 'bytes', d(i).bytes); %#ok<AGROW>
    end
end
end

% ---------------- helper: VDM files that appeared or changed ----------------

function newFiles = diff_vdm_candidates(before, after)
% DIFF_VDM_CANDIDATES  Paths in `after` that are new, or were rewritten.

newFiles = {};
beforePaths = {before.path};

for i = 1:numel(after)
    idx = find(strcmp(after(i).path, beforePaths), 1);
    if isempty(idx)
        newFiles{end+1,1} = after(i).path; %#ok<AGROW>
    elseif after(i).datenum > before(idx).datenum || after(i).bytes ~= before(idx).bytes
        newFiles{end+1,1} = after(i).path; %#ok<AGROW>
    end
end
end

% ---------------- helper: pick newest file ----------------

function f = pick_newest_file(fileList)
assert(~isempty(fileList));
bestT = -Inf;
f = fileList{1};

for i = 1:numel(fileList)
    info = dir(fileList{i});
    if ~isempty(info) && info.datenum > bestT
        bestT = info.datenum;
        f = fileList{i};
    end
end
end
