function join_masks(root, labelsA, labelsB, outLabel)
%JOIN_MASKS  Union of two ROI masks, written as one binary mask per hemisphere.
%
%   join_masks('/Volumes/train_juli/SOLIS/bids')
%   join_masks(someOtherDir)
%   join_masks(someDir, {'lpmErC','rpmErC'}, 'princeErC', 'pmprinceErC')
%
% Walks ROOT recursively, and for every mask it finds of the first label it
% looks for the matching mask of the second label in the same folder, ORs the
% two together and writes the result next to them under OUTLABEL.
%
% The output is BINARY: every voxel in either input, and the overlap, comes
% out as 1. Nothing in it records which of the two regions a voxel came from,
% so the numbers printed per file are the only record of how much each input
% contributed - read them, they are the check that the join did what you meant.
%
% HEMISPHERES STAY SEPARATE. Left is joined with left, right with right, and
% the hemi- entity carries through to the output name. Nothing here produces a
% bilateral mask.
%
% LABEL ARGUMENTS
%   One name  ('princeErC')          used for both hemispheres
%   Two names ({'lpmErC','rpmErC'})  {left, right} - for labels that carry
%                                    the hemisphere in the name, as pmErC does
%
% The defaults are this project's case: lpmErC/rpmErC joined with princeErC
% into pmprinceErC.
%
% NAMING. Only the label is substituted, so everything else in the filename is
% kept, the run- entity included:
%   sub-01s17_ses-01_hemi-left_label-lpmErC_mask.nii
% + sub-01s17_ses-01_hemi-left_label-princeErC_mask.nii
% = sub-01s17_ses-01_hemi-left_label-pmprinceErC_mask.nii
%
% Outputs are overwritten on a re-run: the join is fully determined by its two
% inputs, so there is nothing in an old copy worth keeping.
%
% A file with no partner is reported and skipped. A partner on a different
% voxel grid is reported and skipped too, and the function errors at the end
% so a mismatch cannot pass unnoticed - masks that do not share a grid cannot
% be joined voxel-wise, and that means something is wrong upstream.
%
% Needs only base MATLAB (niftiread/niftiwrite). SPM does not have to be on
% the path.
%
% HEADER. The output carries the first input's grid: same dimensions, same
% voxel sizes, same sform, so it is already aligned with the T2w the two masks
% were drawn on and needs no resampling. niftiwrite writes the sform only, so
% the output has qform_code 0 where the inputs had 2. Nothing is lost - the
% sform is identical and is what SPM, FSL and nibabel read first - but the
% field does differ from the masks the Python extraction wrote.

if nargin < 2 || isempty(labelsA),  labelsA  = {'lpmErC', 'rpmErC'}; end
if nargin < 3 || isempty(labelsB),  labelsB  = 'princeErC';          end
if nargin < 4 || isempty(outLabel), outLabel = 'pmprinceErC';        end

root     = char(root);
outLabel = char(outLabel);
assert(isfolder(root), 'Not a folder: %s', root);

hemis   = {'left', 'right'};
labelsA = per_hemi(labelsA);
labelsB = per_hemi(labelsB);

fprintf('\njoin_masks\n');
fprintf('  root  : %s\n', root);
for h = 1:2
    fprintf('  %-5s : %s + %s -> %s\n', hemis{h}, labelsA{h}, labelsB{h}, outLabel);
end

fprintf('\n%-11s %-8s %-6s %8s %8s %8s %8s\n', ...
        'subject', 'session', 'hemi', 'n_A', 'n_B', 'overlap', 'union');

nWritten  = 0;
nSkipped  = 0;
nMismatch = 0;

for h = 1:2
    A = labelsA{h};
    B = labelsB{h};

    files = find_masks(root, hemis{h}, A);
    for k = 1:numel(files)
        fa  = files{k};
        fb  = strrep(fa, ['label-' A '_mask'], ['label-' B '_mask']);
        fo  = strrep(fa, ['label-' A '_mask'], ['label-' outLabel '_mask']);
        tag = sprintf('%-11s %-8s %-6s', entity(fa, 'sub'), entity(fa, 'ses'), hemis{h});

        if ~isfile(fb)
            fprintf('%s  (no %s mask next to %s - skipped)\n', tag, B, shortname(fa));
            nSkipped = nSkipped + 1;
            continue;
        end

        ia = niftiinfo(fa);
        ib = niftiinfo(fb);
        if ~same_grid(ia, ib)
            fprintf('%s  (GRID MISMATCH with %s - skipped)\n', tag, shortname(fb));
            nMismatch = nMismatch + 1;
            continue;
        end

        a = niftiread(ia) > 0;
        b = niftiread(ib) > 0;
        u = a | b;

        fprintf('%s %8d %8d %8d %8d   %s\n', tag, sum(a(:)), sum(b(:)), ...
                sum(a(:) & b(:)), sum(u(:)), shortname(fo));

        write_mask(u, fo, ia);
        nWritten = nWritten + 1;
    end
end

fprintf('\n  written  : %d\n', nWritten);
fprintf('  skipped  : %d  (no partner mask)\n', nSkipped);
fprintf('  mismatch : %d  (partner on a different grid)\n\n', nMismatch);

if nMismatch > 0
    error('join_masks:gridMismatch', ...
          ['%d mask pair(s) do not share a voxel grid and were not joined. ' ...
           'Both labels have to be in the same space before they can be ORed.'], ...
          nMismatch);
end
end


% ---------------------------------------------------------------------------
function c = per_hemi(x)
% One label name for both hemispheres, or {left, right}.
if ischar(x) || isstring(x)
    c = {char(x), char(x)};
else
    c = cellfun(@char, x(:)', 'UniformOutput', false);
end
assert(numel(c) == 2, 'Give one label name, or {left, right}.');
end


function files = find_masks(root, hemi, label)
% Every <hemi>/<label> mask under root, at any depth. Includes root itself,
% so a flat folder of masks works as well as a BIDS tree.
d = dir(fullfile(root, '**', sprintf('*_hemi-%s_label-%s_mask.nii*', hemi, label)));
d = d(~[d.isdir]);
keep = (endsWith({d.name}, '.nii') | endsWith({d.name}, '.nii.gz')) ...
       & ~startsWith({d.name}, '._');   % macOS resource-fork stubs, not masks
d = d(keep);
files = arrayfun(@(x) fullfile(x.folder, x.name), d, 'UniformOutput', false);
if isempty(files), files = {}; return; end
files = sort(files(:))';
end


function tf = same_grid(ia, ib)
tf = isequal(ia.ImageSize, ib.ImageSize) && ...
     max(abs(ia.Transform.T(:) - ib.Transform.T(:))) < 1e-4;
end


function write_mask(u, outPath, info)
% Write the union as uint8 0/1, on the input's grid, in the input's format.
isGz = endsWith(outPath, '.nii.gz');
if isGz
    base = extractBefore(outPath, strlength(outPath) - 6);
else
    base = extractBefore(outPath, strlength(outPath) - 3);
end

info.Datatype             = 'uint8';
info.BitsPerPixel         = 8;
info.MultiplicativeScaling = 1;   % so a 1 in the file reads back as a 1
info.AdditiveOffset        = 0;
info.DisplayIntensityRange = [0 1];
info.Description           = 'join_masks union';

niftiwrite(uint8(u), char(base), info, 'Compressed', isGz);
end


function v = entity(f, key)
v = regexp(f, ['(?<=[/\\])' key '-[^_/\\]+'], 'match', 'once');
if isempty(v)
    v = regexp(f, [key '-[^_/\\]+'], 'match', 'once');
end
if isempty(v), v = '?'; end
end


function s = shortname(f)
[~, n, e] = fileparts(f);
s = [n e];
end
