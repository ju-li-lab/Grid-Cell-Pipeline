function T = check_roi_coverage(derivRoot, varargin)
% CHECK_ROI_COVERAGE  How much of each ROI mask survives SPM's implicit mask.
%
%   check_roi_coverage(derivRoot)
%   check_roi_coverage(derivRoot, 'Trial', 'sep4_fwhm5', 'CSV', 'roi_coverage.csv')
%
% Answers one question per subject-session-hemisphere: of the voxels in the
% resliced ROI mask, how many would actually carry a value in the GLM, at a
% range of SPM masking thresholds. Run it after preprocessing and before
% committing to a GridCAT run, and run it again after changing a
% preprocessing setting to see whether the change helped.
%
% WHAT IT REPRODUCES
% ------------------
% SPM keeps a voxel only if, in EVERY scan t of the GLM,
%
%       value(voxel, t)  >  mthresh * g_t
%
% where g_t is that scan's global (spm_global: the mean of voxels above
% mean/8). Two details matter and are easy to get wrong:
%
%   - the global is computed PER SCAN, not once. Comparing each voxel's
%     minimum over time against a single global exaggerates the loss,
%     because a globally dimmer volume lowers the voxel and the global
%     together.
%   - the test is across ALL runs that go into the same GLM, not one run.
%
% COLUMNS
% -------
%   vox        voxels in the mask
%   no_data    of those, how many have essentially no signal at all (mean
%              below 1). These are dropout or mask-off-the-brain, and no
%              masking threshold will recover them - worth separating from
%              voxels that are merely dim.
%   mean_sig   mean signal in the mask, raw units
%   ratio_p05  5th percentile of (value / global) over the mask; the
%              threshold you would need to keep 95% of it
%   n@<thr>    voxels surviving at that masking threshold
%
% OPTIONS
%   'Trial'       label written into the CSV, so trials can be compared
%   'Thresholds'  default [0.8 0.7 0.6 0.5 0.4 0.3 0.2]
%   'FuncPrefix'  '' (default) auto-detects: smoothed (su*) if present,
%                 else realigned (u*). Pass one explicitly to override.
%   'CSV'         path to append a row per mask to. Created if absent.
%   'Subjects'    {'sub-01s17', ...} to limit the scan
%
% Needs only base MATLAB (niftiread). SPM does not have to be on the path.

p = inputParser;
addRequired( p, 'derivRoot',  @(x) ischar(x) || isstring(x));
addParameter(p, 'Trial',      '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Thresholds', [0.8 0.7 0.6 0.5 0.4 0.3 0.2], @isnumeric);
addParameter(p, 'FuncPrefix', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'CSV',        '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Subjects',   {}, @iscell);
parse(p, derivRoot, varargin{:});
o = p.Results;

derivRoot = char(derivRoot);
assert(isfolder(derivRoot), 'Not a folder: %s', derivRoot);
TH = sort(o.Thresholds, 'descend');

subs = dir(fullfile(derivRoot, 'sub-*'));
subs = subs([subs.isdir]);
if ~isempty(o.Subjects)
    subs = subs(ismember({subs.name}, o.Subjects));
end
assert(~isempty(subs), 'No sub-* directories in %s', derivRoot);

rows = cell(0, 1);

fprintf('\nROI coverage under SPM''s implicit mask\n');
fprintf('  derivatives : %s\n', derivRoot);
if ~isempty(o.Trial), fprintf('  trial       : %s\n', char(o.Trial)); end
fprintf('\n%-11s %-8s %-6s %5s %8s %9s %9s', ...
        'subject', 'session', 'hemi', 'vox', 'no_data', 'mean_sig', 'ratio_p05');
fprintf('%8.2f', TH);
fprintf('\n');

for a = 1:numel(subs)
    sess = dir(fullfile(derivRoot, subs(a).name, 'ses-*'));
    sess = sess([sess.isdir]);
    for b = 1:numel(sess)
        sdir = fullfile(derivRoot, subs(a).name, sess(b).name);

        [files, usedPrefix] = find_func(sdir, char(o.FuncPrefix));
        if isempty(files)
            fprintf('%-11s %-8s  (no functional files found - skipped)\n', ...
                    subs(a).name, sess(b).name);
            continue;
        end

        [ratioMin, meanImg] = scan_session(files);

        masks = dir(fullfile(sdir, 'anat', 'r*mask.nii'));
        if isempty(masks)
            fprintf('%-11s %-8s  (no resliced ROI masks in anat/ - skipped)\n', ...
                    subs(a).name, sess(b).name);
            continue;
        end

        for k = 1:numel(masks)
            M = niftiread(fullfile(masks(k).folder, masks(k).name));
            if numel(M) ~= numel(ratioMin)
                error(['Mask and functional data are on different grids:\n' ...
                       '  %s  (%d voxels)\n  functional (%d voxels)\n' ...
                       'The mask must be resliced into EPI space - check the coreg stage.'], ...
                       masks(k).name, numel(M), numel(ratioMin));
            end
            idx = find(M(:) > 0.5);
            if isempty(idx), continue; end

            r       = double(ratioMin(idx));
            sig     = double(meanImg(idx));
            nvox    = numel(idx);
            noData  = sum(sig < 1);
            n       = arrayfun(@(th) sum(r > th), TH);
            p05     = prctile(r, 5);

            fprintf('%-11s %-8s %-6s %5d %8d %9.1f %9.2f', ...
                    subs(a).name, sess(b).name, hemi_of(masks(k).name), ...
                    nvox, noData, mean(sig), p05);
            fprintf('%8d', n);
            fprintf('\n');

            rows{end+1, 1} = struct( ...
                'trial',      char(o.Trial), ...
                'subject',    subs(a).name, ...
                'session',    sess(b).name, ...
                'hemi',       hemi_of(masks(k).name), ...
                'prefix',     usedPrefix, ...
                'vox',        nvox, ...
                'no_data',    noData, ...
                'mean_sig',   mean(sig), ...
                'ratio_p05',  p05, ...
                'counts',     n); %#ok<AGROW>
        end
    end
end

T = summarise(rows, TH);

if ~isempty(o.CSV)
    write_csv(char(o.CSV), rows, TH);
end

end


% ---------------------------------------------------------------------------

function [files, prefix] = find_func(sdir, wanted)
% The GLM reads the smoothed data when smoothing ran, the realigned data
% otherwise. Prefer whichever is actually on disk unless told.
if ~isempty(wanted)
    files  = dir(fullfile(sdir, 'func', [wanted '*_bold.nii']));
    prefix = wanted;
    return;
end
for prefix = {'su', 'u'}
    files = dir(fullfile(sdir, 'func', [prefix{1} '*_bold.nii']));
    files = files(~startsWith({files.name}, 'mean'));   % not the mean image
    if ~isempty(files), prefix = prefix{1}; return; end
end
files  = [];
prefix = '';
end


function [ratioMin, meanImg] = scan_session(files)
% Per voxel: the smallest value/global ratio seen in any scan of any run,
% and the mean signal across everything. One volume at a time, so a long
% session does not need the whole 4D series in memory twice.
ratioMin = [];
total    = [];
nScans   = 0;

for k = 1:numel(files)
    V  = niftiread(fullfile(files(k).folder, files(k).name));
    nT = size(V, 4);
    Vr = reshape(single(V), [], nT);
    clear V

    if isempty(ratioMin)
        ratioMin = inf(size(Vr, 1), 1, 'single');
        total    = zeros(size(Vr, 1), 1, 'double');
    elseif size(Vr, 1) ~= numel(ratioMin)
        error('Run %s has a different voxel grid from the first run in this session.', ...
              files(k).name);
    end

    for t = 1:nT
        v = Vr(:, t);
        m = mean(v);
        g = mean(v(v > m/8));          % spm_global
        if ~isfinite(g) || g <= 0, continue; end
        ratioMin = min(ratioMin, v / g);
        total    = total + double(v);
        nScans   = nScans + 1;
    end
    clear Vr
end

meanImg = total / max(nScans, 1);
end


function h = hemi_of(name)
lower_name = lower(name);
if     contains(lower_name, 'left'),  h = 'left';
elseif contains(lower_name, 'right'), h = 'right';
elseif contains(lower_name, 'bilat'), h = 'bilat';
else,                                 h = 'other';
end
end


function T = summarise(rows, TH)
if isempty(rows)
    fprintf('\nNothing measured.\n');
    T = table();
    return;
end
C   = cell2mat(cellfun(@(r) r.counts(:)',  rows, 'UniformOutput', false));
vox = cell2mat(cellfun(@(r) r.vox,         rows, 'UniformOutput', false));
pct = 100 * C ./ vox;

fprintf('\n%d masks measured\n\n', numel(rows));
fprintf('%-9s %9s %8s %18s\n', 'mthresh', 'median %', 'worst %', 'masks losing >half');
for i = 1:numel(TH)
    fprintf('%-9.2f %8.0f%% %7.0f%% %12d / %-4d\n', ...
            TH(i), median(pct(:,i)), min(pct(:,i)), sum(pct(:,i) < 50), size(pct,1));
end
fprintf(['\nPick the lowest threshold that is still on the flat part of this\n' ...
         'curve. If a mask has a high no_data count, no threshold will fix it -\n' ...
         'that is dropout or a misplaced mask, not conservative masking.\n']);

T = table(TH(:), median(pct,1)', min(pct,[],1)', sum(pct<50,1)', ...
          'VariableNames', {'mthresh','median_pct','worst_pct','masks_below_half'});
end


function write_csv(path, rows, TH)
newFile = ~isfile(path);
fid = fopen(path, 'a');
assert(fid > 0, 'Cannot write %s', path);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

if newFile
    fprintf(fid, 'trial,subject,session,hemi,prefix,vox,no_data,mean_sig,ratio_p05');
    fprintf(fid, ',n_%.2f', TH);
    fprintf(fid, ',pct_%.2f', TH);
    fprintf(fid, '\n');
end
for i = 1:numel(rows)
    r = rows{i};
    fprintf(fid, '%s,%s,%s,%s,%s,%d,%d,%.2f,%.4f', ...
            r.trial, r.subject, r.session, r.hemi, r.prefix, ...
            r.vox, r.no_data, r.mean_sig, r.ratio_p05);
    fprintf(fid, ',%d', r.counts);
    fprintf(fid, ',%.1f', 100 * r.counts / r.vox);
    fprintf(fid, '\n');
end
fprintf('\nAppended %d rows to %s\n', numel(rows), path);
end
