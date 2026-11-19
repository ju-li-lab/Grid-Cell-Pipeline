function T = compare_roi_sessions(derivRoot, refSes, movSes, subjects)
%COMPARE_ROI_SESSIONS  Put one session's ROI masks on another session's EPI
%grid and measure how much they overlap and how far they moved.
%
%   compare_roi_sessions('/Volumes/train_juli/spm-preproc-unwarp1-stc2nd')
%   compare_roi_sessions(derivRoot, 'ses-01', 'ses-03')
%   compare_roi_sessions(derivRoot, 'ses-01', 'ses-02', {'sub-05s17'})
%
% Per subject, the movSes session mean is coregistered to the refSes session
% mean (estimate only) with the movSes masks riding along, and the masks are
% then resliced once onto the refSes EPI grid, the grid GridCAT reads. They
% start from the native-resolution masks in movSes/anat, whose headers the
% coreg stage already moved onto the movSes mean, and are resliced with the
% options from refSes's spm_options.m snapshot, so they are made exactly the
% way the pipeline made the refSes r*_mask.nii they are compared against.
%
% Columns:
%   n_ref      voxels in the refSes mask
%   n_mov_own  voxels in movSes's own resliced mask (what GridCAT got there)
%   n_mov      voxels in the movSes mask on the refSes grid
%   n_both     voxels in both
%   dice       2 * n_both / (n_ref + n_mov)
%   dx dy dz   movSes centroid minus refSes centroid, mm, along the refSes
%              scanner axes (x left->right, y back->front, z down->up)
%
% Writes only into <derivRoot>/roi_<movSes>_to_<refSes>/:
%   <sub>/r<mean>          movSes mean on the refSes grid; lay it over the
%                          refSes mean first to check the coregistration
%   <sub>/r<mask>          movSes masks on the refSes grid
%   <sub>/overlap_<mask>   1 = refSes only, 2 = movSes only, 3 = both
%   roi_<movSes>_to_<refSes>.csv
%
% Needs SPM on the path.

if nargin < 2 || isempty(refSes), refSes = 'ses-01'; end
if nargin < 3 || isempty(movSes), movSes = 'ses-02'; end
if nargin < 4 || isempty(subjects)
    d = dir(fullfile(derivRoot, 'sub-*'));
    subjects = {d([d.isdir]).name};
end

derivRoot = char(derivRoot);
outRoot   = fullfile(derivRoot, sprintf('roi_%s_to_%s', movSes, refSes));

spm('defaults', 'fmri');
spm_jobman('initcfg');

rows = {};
for s = 1:numel(subjects)
    sub    = subjects{s};
    refDir = fullfile(derivRoot, sub, refSes);
    movDir = fullfile(derivRoot, sub, movSes);

    refMean  = list(fullfile(refDir, 'func', 'mean*.nii'));
    movMean  = list(fullfile(movDir, 'func', 'mean*.nii'));
    movMasks = list(fullfile(movDir, 'anat', [sub '_' movSes '*_mask.nii']));
    if isempty(refMean) || isempty(movMean) || isempty(movMasks)
        fprintf('%-11s (no session mean or no masks in %s/%s - skipped)\n', sub, refSes, movSes);
        continue;
    end
    refMean = refMean{1};

    % the refSes mask each movSes mask is compared against, matched on hemi and label
    refMasks = cell(size(movMasks));
    for k = 1:numel(movMasks)
        suffix = regexp(movMasks{k}, '_hemi-[^/\\]*$', 'match', 'once');
        hit    = list(fullfile(refDir, 'anat', ['r' sub '_' refSes '*' suffix]));
        if ~isempty(hit), refMasks{k} = hit{1}; end
    end
    keep     = ~cellfun(@isempty, refMasks);
    movMasks = movMasks(keep);
    refMasks = refMasks(keep);
    if isempty(movMasks)
        fprintf('%-11s (no %s mask matches a %s mask - skipped)\n', sub, movSes, refSes);
        continue;
    end

    % SPM rewrites the headers it estimates, so it works on copies
    outDir = fullfile(outRoot, sub);
    if ~isfolder(outDir), mkdir(outDir); end
    meanCopy   = copy_into(movMean{1}, outDir);
    maskCopies = cellfun(@(f) copy_into(f, outDir), movMasks, 'UniformOutput', false);

    addpath(refDir, '-begin');      % the snapshot that made the refSes masks
    OPT = spm_options();
    rmpath(refDir);

    mb = {};
    mb{1}.spm.spatial.coreg.estimate.ref      = {refMean};
    mb{1}.spm.spatial.coreg.estimate.source   = {meanCopy};
    mb{1}.spm.spatial.coreg.estimate.other    = maskCopies(:);
    mb{1}.spm.spatial.coreg.estimate.eoptions = OPT.coreg.eoptions;
    spm_jobman('run', mb);

    reslice(refMean, maskCopies, OPT.coreg.roptions);
    ro        = OPT.coreg.roptions;
    ro.interp = 4;                   % the mean is not a label image
    reslice(refMean, {meanCopy}, ro);

    for k = 1:numel(maskCopies)
        Vref = spm_vol(refMasks{k});
        A    = spm_read_vols(Vref) > 0.5;
        B    = spm_read_vols(spm_vol(spm_file(maskCopies{k}, 'prefix', 'r'))) > 0.5;
        own  = spm_read_vols(spm_vol(spm_file(movMasks{k}, 'prefix', 'r'))) > 0.5;

        d    = centroid(B, Vref.mat) - centroid(A, Vref.mat);
        hemi = regexp(maskCopies{k}, '(?<=_hemi-)[^_]+', 'match', 'once');
        lab  = regexp(maskCopies{k}, '(?<=_label-)[^_]+', 'match', 'once');
        dice = 2 * nnz(A & B) / (nnz(A) + nnz(B));

        rows(end+1, :) = {sub, hemi, lab, nnz(A), nnz(own), nnz(B), nnz(A & B), ...
                          dice, d(1), d(2), d(3), norm(d)}; %#ok<AGROW>

        V         = rmfield(Vref, 'private');
        V.fname   = fullfile(outDir, ['overlap_' spm_file(maskCopies{k}, 'filename')]);
        V.dt      = [spm_type('uint8') spm_platform('bigend')];
        V.pinfo   = [1; 0; 0];
        V.descrip = sprintf('1=%s 2=%s 3=both', refSes, movSes);
        spm_write_vol(V, double(A) + 2 * double(B));
    end
end

T = cell2table(rows, 'VariableNames', {'subject', 'hemi', 'label', 'n_ref', ...
    'n_mov_own', 'n_mov', 'n_both', 'dice', 'dx', 'dy', 'dz', 'shift_mm'});
if ~isempty(rows)
    fprintf('\n%s masks on the %s EPI grid\n\n', movSes, refSes);
    fprintf('%-11s %-6s %6s %9s %6s %7s %6s %6s %6s %6s %8s\n', 'subject', 'hemi', ...
            'n_ref', 'n_mov_own', 'n_mov', 'n_both', 'dice', 'dx', 'dy', 'dz', 'shift');
    for r = 1:size(rows, 1)
        fprintf('%-11s %-6s %6d %9d %6d %7d %6.2f %6.1f %6.1f %6.1f %5.1f mm\n', rows{r, [1 2 4:12]});
    end
    csv = fullfile(outRoot, sprintf('roi_%s_to_%s.csv', movSes, refSes));
    writetable(T, csv);
    fprintf('\nWrote %s\n', csv);
end
end


function f = list(pattern)
d = dir(pattern);
d = d(~[d.isdir] & ~startsWith({d.name}, '._'));   % macOS stubs on the external drive
f = arrayfun(@(x) fullfile(x.folder, x.name), d, 'UniformOutput', false);
end


function out = copy_into(src, destDir)
out = fullfile(destDir, spm_file(src, 'filename'));
copyfile(src, out, 'f');
fileattrib(out, '+w');
end


function reslice(ref, src, roptions)
roptions.prefix = 'r';
mb = {};
mb{1}.spm.spatial.coreg.write.ref      = {ref};
mb{1}.spm.spatial.coreg.write.source   = src(:);
mb{1}.spm.spatial.coreg.write.roptions = roptions;
spm_jobman('run', mb);
end


function c = centroid(M, mat)
[i, j, k] = ind2sub(size(M), find(M));
c = mat * [mean([i j k], 1)'; 1];
c = c(1:3)';
end
