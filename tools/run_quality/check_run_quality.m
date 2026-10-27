function T = check_run_quality(derivRoot, varargin)
% CHECK_RUN_QUALITY  Per-run data quality, with each metric's reliability measured.
%
%   check_run_quality(derivRoot)
%   check_run_quality(derivRoot, 'Trial','sep4_fwhm5', 'CSV','run_quality.csv')
%
% Which runs are worth keeping? This reports four things per run and, at the
% end, the SPLIT-HALF RELIABILITY of each one measured on your own data - so
% you can see which numbers you are entitled to act on.
%
% WHY THE RELIABILITY CHECK IS NOT OPTIONAL
% -----------------------------------------
% A per-run number is only useful if the same run scores the same twice. The
% task-fit R2 below was proposed as a quality metric and then failed exactly
% that test on this dataset: split-half r = 0.10, i.e. the first and second
% half of one run disagree almost completely. Acting on it would mean
% dropping runs at random. It is still computed and reported, because on a
% dataset with a stronger task response it may well be reliable - but the
% tool measures that for you rather than assuming it.
%
% THE METRICS
%   tsnr_brain   median temporal SNR over the brain. Model-free, and normally
%                the most reliable number here.
%   tsnr_roi     the same inside the ROI mask. More relevant than tsnr_brain
%                for an ROI analysis: the ErC can be poor while the brain
%                looks fine.
%   mean_fd      framewise displacement (Power), mm. |dtrans| + 50*|drot|.
%   pct_fd_hi    share of volumes above FdThreshold.
%   dvars        median standardised temporal derivative of the signal.
%                Catches spikes and dropout that motion parameters miss.
%   task_r2      partial R2 of the task regressors (rotate/stand/translate)
%                in the 500 voxels that respond most, over and above motion
%                and drift. See the warning above.
%
% TWO THINGS THAT MAKE task_r2 HONEST, AND MATTER
%   - voxels are chosen LEAVE-ONE-RUN-OUT, from a contrast computed on the
%     OTHER runs of that session. Choosing them from all runs inflated this
%     by +131% on the test dataset, and by 23x for individual runs.
%   - motion regressors are in the nuisance model. Leaving them out inflated
%     it by a further +62%, because a run whose motion follows the task
%     otherwise scores well for the wrong reason.
%
% OPTIONS
%   'GLMRoot'     GLM output holding <sub>_<ses>_GLM1 (needed for task_r2;
%                 without it the model-free metrics are still reported)
%   'EventDir'    EventData directory
%   'ROIDir'      ROI masks, for tsnr_roi
%   'Sessions'    {} = all. Unlike the grid-consistency tool this defaults to
%                 every session: data quality is not affected by an
%                 intervention the way grid signal might be.
%   'Subjects'    {} = all
%   'FuncPrefix'  '' auto-detects su* then u*
%   'FdThreshold' 0.5 mm
%   'NTop'        500 voxels for task_r2
%   'Trial'       label for the CSV
%   'CSV'         path to append to
%
% Needs only base MATLAB. SPM and GridCAT are not required.

p = inputParser;
addRequired( p,'derivRoot',  @(x)ischar(x)||isstring(x));
addParameter(p,'GLMRoot',    '', @(x)ischar(x)||isstring(x));
addParameter(p,'EventDir',   '', @(x)ischar(x)||isstring(x));
addParameter(p,'ROIDir',     '', @(x)ischar(x)||isstring(x));
addParameter(p,'Sessions',   {}, @iscell);
addParameter(p,'Subjects',   {}, @iscell);
addParameter(p,'FuncPrefix', '', @(x)ischar(x)||isstring(x));
addParameter(p,'FdThreshold',0.5, @isnumeric);
addParameter(p,'NTop',       500, @isnumeric);
addParameter(p,'Trial',      '', @(x)ischar(x)||isstring(x));
addParameter(p,'CSV',        '', @(x)ischar(x)||isstring(x));
parse(p, derivRoot, varargin{:}); o = p.Results;

derivRoot = char(derivRoot);
assert(isfolder(derivRoot), 'Not a folder: %s', derivRoot);
haveGLM = ~isempty(o.GLMRoot) && isfolder(char(o.GLMRoot)) && ~isempty(o.EventDir);

fprintf('\nPer-run data quality\n  derivatives : %s\n', derivRoot);
if haveGLM, fprintf('  GLM1 output : %s\n', char(o.GLMRoot));
else,       fprintf('  (no GLMRoot/EventDir given - task_r2 skipped)\n'); end
fprintf('\n%-11s %-8s %-7s %8s %8s %8s %8s %8s %9s\n', ...
        'subject','session','run','tsnr_br','tsnr_roi','mean_fd','pct_fd','dvars','task_r2');

rows = {};
subs = dir(fullfile(derivRoot,'sub-*')); subs = subs([subs.isdir]);
if ~isempty(o.Subjects), subs = subs(ismember({subs.name}, o.Subjects)); end

for a = 1:numel(subs)
  sess = dir(fullfile(derivRoot, subs(a).name, 'ses-*')); sess = sess([sess.isdir]);
  if ~isempty(o.Sessions), sess = sess(ismember({sess.name}, o.Sessions)); end
  for b = 1:numel(sess)
    sub = subs(a).name; ses = sess(b).name;
    fdir = fullfile(derivRoot, sub, ses, 'func');
    [files, pre] = find_runs(fdir, char(o.FuncPrefix));
    if isempty(files), continue; end

    roi = [];
    if ~isempty(o.ROIDir)
      m = dir(fullfile(char(o.ROIDir), sprintf('r%s_%s_hemi-bilat*_mask.nii', sub, ses)));
      if ~isempty(m), roi = niftiread(fullfile(m(1).folder, m(1).name)) > 0.5; end
    end

    con = {};
    if haveGLM
      con = loo_contrasts(fullfile(char(o.GLMRoot), sprintf('%s_%s_GLM1', sub, ses)), numel(files));
    end

    for r = 1:numel(files)
      Y4 = double(niftiread(fullfile(files(r).folder, files(r).name)));
      sz = size(Y4); nv = sz(4);
      Y  = reshape(Y4, [], nv)'; clear Y4

      mu = mean(Y,1); sd = std(Y,0,1);
      brain = mu > 0.2*median(mu(mu>0));            % crude but adequate brain mask
      tsnr_br  = median(mu(brain)./max(sd(brain),eps));
      tsnr_roi = NaN;
      if ~isempty(roi) && numel(roi) == size(Y,2)
          k = find(roi(:)); tsnr_roi = median(mu(k)./max(sd(k),eps));
      end

      dv = sqrt(mean(diff(Y(:,brain)).^2, 2));
      dvars = median(dv) / median(mu(brain)) * 100;  % percent of mean signal

      [mfd, pfd, M] = motion_stats(fdir, files(r).name, pre, nv, o.FdThreshold);

      r2 = NaN; r2a = NaN; r2b = NaN;
      if haveGLM && numel(con) >= r && ~isempty(con{r}) && ~isempty(M)
        task = task_regressors(char(o.EventDir), sub, ses, r, nv, run_tr(files(r)));
        if ~isempty(task)
          v = con{r}(:); v(~isfinite(v)) = -inf;
          [~,ord] = sort(v,'descend'); idx = ord(1:min(o.NTop, numel(ord)));
          y = mean(Y(:,idx), 2);
          TR = run_tr(files(r)); h1 = 1:floor(nv/2); h2 = floor(nv/2)+1:nv;
          r2  = partial_r2(y,      task,         [dctmat(nv,TR,128)         M        ]);
          r2a = partial_r2(y(h1),  task(h1,:),   [dctmat(numel(h1),TR,128)  M(h1,:)  ]);
          r2b = partial_r2(y(h2),  task(h2,:),   [dctmat(numel(h2),TR,128)  M(h2,:)  ]);
        end
      end
      clear Y

      fprintf('%-11s %-8s task-%-2d %8.1f %8.1f %8.3f %7.0f%% %8.2f %9.4f\n', ...
              sub, ses, r, tsnr_br, tsnr_roi, mfd, 100*pfd, dvars, r2);
      rows(end+1,:) = {char(o.Trial),sub,ses,r,tsnr_br,tsnr_roi,mfd,pfd,dvars,r2,r2a,r2b}; %#ok<AGROW>
    end
  end
end

T = summarise(rows);
if ~isempty(o.CSV), write_csv(char(o.CSV), rows); end
end


% ---------------------------------------------------------------------------

function [files, prefix] = find_runs(fdir, wanted)
if ~isempty(wanted)
    files = dir(fullfile(fdir, [wanted '*_bold.nii'])); prefix = wanted; return;
end
for pre = {'su','u'}
    files = dir(fullfile(fdir, [pre{1} 'sub*_bold.nii']));
    if ~isempty(files), prefix = pre{1}; return; end
end
files = []; prefix = '';
end

function TR = run_tr(f)
j = regexprep(fullfile(f.folder,f.name), '^(.*/)su?(sub-.*)\.nii$', '$1$2.json');
TR = 2;
if isfile(j)
    d = jsondecode(fileread(j));
    if isfield(d,'RepetitionTime'), TR = d.RepetitionTime; end
end
end

function con = loo_contrasts(gdir, nruns)
% translate - stand, computed on the OTHER runs. Leave-one-run-out, so a run
% never helps choose the voxels it is then scored in.
con = {};
if ~isfile(fullfile(gdir,'SPM.mat')), return; end
L = load(fullfile(gdir,'SPM.mat')); nm = L.SPM.xX.name;
B = struct('t',{{}},'s',{{}});
for s = 1:nruns
    kt = find(strcmp(nm, sprintf('Sn(%d) GridEvent-translate*bf(1)', s)));
    ks = find(strcmp(nm, sprintf('Sn(%d) stand*bf(1)', s)));
    if isempty(kt) || isempty(ks), return; end
    B.t{s} = double(niftiread(fullfile(gdir, sprintf('beta_%04d.nii', kt))));
    B.s{s} = double(niftiread(fullfile(gdir, sprintf('beta_%04d.nii', ks))));
end
for r = 1:nruns
    other = setdiff(1:nruns, r); c = zeros(size(B.t{1}));
    for s = other, c = c + (B.t{s} - B.s{s})/numel(other); end
    con{r} = c; %#ok<AGROW>
end
end

function [mfd, pfd, M] = motion_stats(fdir, fname, pre, nv, thr)
stem = regexprep(fname, ['^' pre '(sub-.*)_bold\.nii$'], '$1');
f = dir(fullfile(fdir, ['rp_' stem '_bold.txt']));
if isempty(f), f = dir(fullfile(fdir, 'rp_*.txt')); end
mfd = NaN; pfd = NaN; M = [];
if isempty(f), return; end
M = readmatrix(fullfile(f(1).folder, f(1).name));
if size(M,1) < nv || size(M,2) < 6, M = []; return; end
M = M(1:nv, 1:6);
fd = [0; sum(abs(diff(M(:,1:3))),2) + 50*sum(abs(diff(M(:,4:6))),2)];
mfd = mean(fd); pfd = mean(fd > thr);
end

function D = task_regressors(evDir, sub, ses, run, nv, TR)
D = [];
f = fullfile(evDir, sprintf('%s_%s_task-%d_EventData.txt', sub, ses, run));
if ~isfile(f), return; end
E = readmatrix(f, 'OutputType','string', 'Delimiter',',');
typ = E(:,1); ons = str2double(E(:,2)); dur = max(str2double(E(:,3)), 0.1);
dt = 0.1; th = 0:dt:32; g = @(x,k) x.^(k-1).*exp(-x)./gamma(k);
hrf = g(th,6) - g(th,16)/6; hrf = hrf/sum(hrf);
types = unique(typ, 'stable'); D = zeros(nv, numel(types));
for q = 1:numel(types)
    sel = typ == types(q); oo = ons(sel); uu = dur(sel);
    hi = zeros(round(nv*TR/dt) + 800, 1);
    for e = 1:numel(oo)
        a = round(oo(e)/dt)+1; b = a + round(uu(e)/dt) - 1;
        a = max(a,1); b = min(b, numel(hi));
        if b >= a, hi(a:b) = 1; end
    end
    cv = conv(hi, hrf);
    D(:,q) = cv(round(((0:nv-1)*TR + TR/2)/dt) + 1);
end
end

function r2 = partial_r2(y, D, N)
y = y(:); y = y - N*(N\y); Dn = D - N*(N\D);
b = Dn\y; r2 = 1 - var(y - Dn*b)/var(y);
end

function K = dctmat(N, TR, cut)
n = max(fix(2*(N*TR)/cut + 1), 1); K = zeros(N, n);
for j = 1:n, K(:,j) = cos(pi*(0:N-1)'*(j-1)/(N-1)); end
end

function T = summarise(rows)
if isempty(rows), fprintf('\nNothing measured.\n'); T = table(); return; end
nm = {'tsnr_brain','tsnr_roi','mean_fd','dvars','task_r2'};
col = [5 6 7 9 10];
fprintf('\n%d runs measured.\n\n', size(rows,1));
fprintf('%-12s %10s %10s %10s\n', 'metric', 'median', 'min', 'max');
for i = 1:numel(nm)
    v = cell2mat(rows(:,col(i)));
    fprintf('%-12s %10.3f %10.3f %10.3f\n', nm{i}, median(v,'omitnan'), min(v), max(v));
end

a = cell2mat(rows(:,11)); b = cell2mat(rows(:,12));
ok = isfinite(a) & isfinite(b);
fprintf('\n---- can you act on task_r2? ----\n');
n = sum(ok);
if n < 10
    % A correlation on a handful of runs is itself unreliable: on one subset
    % of 6 runs this read r = 0.89 where the full 18 gave r = 0.10. Refuse to
    % pronounce rather than invite a decision from noise.
    fprintf('  only %d runs with a task fit - too few to judge reliability.\n', n);
    fprintf('  Run over more subjects/sessions before trusting task_r2 either way.\n');
else
    r  = corr(a(ok), b(ok));
    sb = 2*r/(1+r);
    % Fisher-z 95%% interval, so the uncertainty in r itself is visible
    z  = atanh(max(min(r,0.999),-0.999)); se = 1/sqrt(n-3);
    lo = tanh(z-1.96*se); hi = tanh(z+1.96*se);
    fprintf('  split-half correlation : r = %.3f   95%% CI [%.2f %.2f]   (%d runs)\n', r, lo, hi, n);
    fprintf('  Spearman-Brown         : r = %.3f\n', sb);
    sbhi = 2*hi/(1+hi);
    if sbhi < 0.7
        fprintf(['  -> NOT RELIABLE on this dataset. The two halves of one run\n' ...
                 '     disagree, so a run''s score is mostly noise. Do not drop\n' ...
                 '     runs on it. Use tsnr_roi and mean_fd instead.\n']);
    elseif sb < 0.7
        fprintf(['  -> NO EVIDENCE it works. The point estimate is poor; the\n' ...
                 '     interval is too wide to rule out moderate reliability,\n' ...
                 '     but you have no grounds to drop runs on it either way.\n' ...
                 '     Use tsnr_roi and mean_fd, which need no such check.\n']);
    else
        fprintf('  -> reliable enough to act on.\n');
    end
end
fprintf(['\n  tsnr_* and mean_fd are model-free and normally reliable.\n' ...
         '  tsnr_roi is the one that matters for an ROI analysis: the brain\n' ...
         '  can look fine while the ErC does not.\n']);

T = cell2table(rows(:,2:10), 'VariableNames', ...
    {'subject','session','run','tsnr_brain','tsnr_roi','mean_fd','pct_fd_hi','dvars','task_r2'});
end

function write_csv(path, rows)
newFile = ~isfile(path); fid = fopen(path,'a');
assert(fid > 0, 'Cannot write %s', path);
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
if newFile
    fprintf(fid, ['trial,subject,session,run,tsnr_brain,tsnr_roi,mean_fd,' ...
                  'pct_fd_hi,dvars,task_r2,task_r2_half1,task_r2_half2\n']);
end
for i = 1:size(rows,1)
    fprintf(fid, '%s,%s,%s,%d,%.3f,%.3f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f\n', rows{i,:});
end
fprintf('\nAppended %d rows to %s\n', size(rows,1), path);
end
