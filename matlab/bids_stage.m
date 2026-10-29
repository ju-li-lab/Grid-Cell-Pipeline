function out = bids_stage(srcFile, destDir, newName)
% BIDS_STAGE  Copy one input file into the derivatives, ready for SPM.
%
%   out = bids_stage(srcFile, destDir)
%   out = bids_stage(srcFile, destDir, 'sub-01_ses-01_task-run1_bold.nii')
%
%   - .nii.gz is decompressed on the way in, so the raw BIDS directory never
%     gains an uncompressed twin (which is what ensure_nii used to do).
%   - The JSON sidecar comes along, so the validation and provenance steps can
%     read the acquisition parameters from the derivative alone.
%   - The .nii is ALWAYS re-copied, never reused. See below.
%
%   Returns the path of the staged .nii.
%
%   Copying is not an optimisation the pipeline can skip: SPM's Realign &
%   Unwarp and Coregister write the estimated transforms into the *input*
%   image headers. A symlink would write straight through to the raw BIDS file.
%
%   That same fact is why a staged .nii is never reused. This function used to
%   skip the copy when the staged file was no older than its source, which
%   looks safe and is not: SPM rewriting the header updates the staged file's
%   modification time, so from the first run onwards it is permanently newer
%   than the raw file and the check says "up to date" forever. Whatever
%   transform SPM last wrote then becomes the starting point of the next run.
%
%   One session failed exactly this way. A coregistration diverged and put the
%   T2w 145 mm from the EPI; every later run reused that header, so the
%   coregistration could never recover and the ROI masks resliced to zero
%   voxels no matter what else was changed. The mtime of a file SPM has
%   written to cannot tell you whether it is still a faithful copy.
%
%   The JSON sidecar is a different matter and is still skipped when current:
%   nothing in the pipeline writes to it, so its timestamp means what it says.
%
%   See also BIDS_DERIV, BIDS_RESET_DERIV.

if isempty(srcFile)
    out = '';
    return;
end

assert(isfile(srcFile), 'bids_stage:MissingSource', 'Cannot stage a file that does not exist:\n  %s', srcFile);

if ~isfolder(destDir)
    [ok, msg] = mkdir(destDir);
    assert(ok, 'bids_stage:CannotCreate', 'Could not create %s\n  %s', destDir, msg);
end

isGz = endsWith(lower(srcFile), '.nii.gz');

if nargin < 3 || isempty(newName)
    [~, base, ext] = fileparts(srcFile);
    if isGz
        newName = base;              % '...nii' — the .gz is dropped by decompressing
    else
        newName = [base ext];
    end
end
newName = regexprep(newName, '\.gz$', '');

out = fullfile(destDir, newName);

if isGz
    fprintf('  [stage] decompress %s -> %s\n', shortname(srcFile), newName);
    tmpDir = tempname;
    mkdir(tmpDir);
    cleanupTmp = onCleanup(@() rmdir_safe(tmpDir));
    produced = gunzip(srcFile, tmpDir);
    assert(~isempty(produced), 'bids_stage:GunzipFailed', 'Could not decompress %s', srcFile);
    if isfile(out), delete(out); end
    movefile(produced{1}, out);
else
    fprintf('  [stage] copy %s -> %s\n', shortname(srcFile), newName);
    if isfile(out), delete(out); end
    copyfile(srcFile, out);
end

% Copied files inherit the source's read-only bit on some filesystems; SPM has
% to be able to rewrite the header it just estimated.
try
    fileattrib(out, '+w');
catch
end

stage_sidecar(srcFile, destDir, newName);
end

% ============================== helpers ==============================

function stage_sidecar(srcFile, destDir, newName)
% Bring the JSON sidecar across under the staged file's name.
srcJson = regexprep(srcFile, '\.nii(\.gz)?$', '.json');
if ~isfile(srcJson), return; end

destJson = fullfile(destDir, regexprep(newName, '\.nii(\.gz)?$', '.json'));
if is_current(destJson, srcJson), return; end

copyfile(srcJson, destJson);
try
    fileattrib(destJson, '+w');
catch
end
end

function tf = is_current(dest, src)
% Only ever asked about JSON sidecars. Do not reach for this to skip copying
% an image: SPM writes to staged images, which makes their timestamps lie.
tf = false;
if ~isfile(dest), return; end
dDest = dir(dest);
dSrc  = dir(src);
if isempty(dDest) || isempty(dSrc), return; end
tf = dDest.datenum >= dSrc.datenum;
end

function rmdir_safe(d)
if isfolder(d)
    try
        rmdir(d, 's');
    catch
    end
end
end

function s = shortname(p)
[~, b, e] = fileparts(p);
s = [b e];
end
