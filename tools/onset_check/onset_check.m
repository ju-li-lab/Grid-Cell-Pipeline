function T = onset_check(derivRoot, eventDir, varargin)
% ONSET_CHECK  Does event time 0 land on the first saved volume?
%
%   onset_check(derivRoot, eventDir)
%   onset_check(derivRoot, eventDir, 'PerSlice', true, 'SliceTiming', st)
%
% Answers one question and nothing else: by how many seconds are the event
% onsets offset from the functional series? It does NOT look for a grid code,
% and it must not - the grid signal is the thing under investigation, so a
% timing test that depended on it could not tell "the clock is wrong" apart
% from "there is no grid". This uses the task's visual structure instead,
% which is there whether or not the entorhinal cortex does anything.
%
% HOW IT WORKS
% ------------
% The rotate / stand / translate events tile the movement phase continuously
% with durations from 0.3 s to 6.4 s, so a design built from them carries
% strong, irregularly spaced structure in the band the haemodynamic response
% passes. Build that design, convolve it with the canonical HRF, slide it
% against the data in steps of one microtime bin, and record the F statistic
% of the task block at every shift. The curve peaks at the true offset.
%
% Voxels are picked by the standard deviation of their nuisance-residualised
% time course. That criterion uses no timing information whatsoever, so it
% cannot pull the peak toward any particular shift. In a passive-navigation
% task the voxels it picks are overwhelmingly visual.
%
% READING THE RESULT  -  the sign matters
% ---------------------------------------
% delta_hat is where the curve peaks, in seconds. The correction is always
%
%       corrected_onset = onset - delta_hat
%
%   delta_hat ~ 0        event time 0 IS the first saved volume. Nothing to do.
%   delta_hat ~ +k*TR    the trigger fired k volumes before the first saved
%                        volume, i.e. k dummies were discarded after the
%                        trigger. Subtract k*TR from every onset.
%   delta_hat ~ -k*TR    the onsets sit k volumes late. Add k*TR.
%
% Whole-volume offsets are what this resolves, and it resolves them
% decisively: per-run precision is a few hundred milliseconds, and since the
% offset is one constant shared by every run the pooled estimate is far
% tighter than one TR. It CANNOT settle sub-second questions such as
% microtimeOnset 8 against 9, because the regional haemodynamic response
% departs from the canonical shape by more than that. Read anything within
% +/- 1 s of zero as zero.
%
% SLICE TIMING
% ------------
% Default mode reads the slice-time-corrected series and samples the design
% at MicrotimeOnset of MicrotimeRes. With MicrotimeOnset = 9, MicrotimeRes =
% 16 and slicetime.ref_fraction = 0.5 those two reference points coincide
% exactly, so delta_hat is the event-clock offset on its own.
%
% 'PerSlice' instead reads a series with NO slice timing correction and fits
% one curve per slice. Each slice's own acquisition time then appears in its
% delta_hat, so regressing delta_hat on SliceTiming should return slope 1,
% and the event-clock offset is  intercept + (MicrotimeOnset-1)*TR/MicrotimeRes.
% That is the independent check: it confirms the offset, the SliceTiming
% vector and the reference fraction in a single fit, assuming none of them.

p = inputParser;
p.addParameter('Subjects',       {});
p.addParameter('Sessions',       {});
p.addParameter('Tasks',          [1 2 3]);
p.addParameter('Prefix',         'sua');
p.addParameter('RawPrefix',      '');
p.addParameter('TR',             2.06);
p.addParameter('MicrotimeRes',   16);
p.addParameter('MicrotimeOnset', 9);
p.addParameter('HPFSec',         128);
p.addParameter('RangeTR',        5);
p.addParameter('KeepPct',        2);
p.addParameter('PerSlice',       false);
p.addParameter('SliceTiming',    []);   % vector of seconds, or a BOLD JSON path
p.addParameter('OutDir',         '');
p.addParameter('CSV',            '');
p.parse(varargin{:});
o = p.Results;

TR = o.TR; Tb = o.MicrotimeRes; T0 = o.MicrotimeOnset; dt = TR/Tb;
ref = (T0 - 1) * dt;
kMax = round(o.RangeTR * TR / dt);
kRange = -kMax:kMax;
prefix = o.Prefix;
if o.PerSlice, prefix = o.RawPrefix; end

outDir = o.OutDir;
if isempty(outDir), outDir = fullfile(fileparts(derivRoot), 'onset_check'); end
if ~isfolder(outDir), mkdir(outDir); end
csvPath = o.CSV;
if isempty(csvPath), csvPath = fullfile(outDir, 'onset_check.csv'); end

fprintf('\nonset_check\n');
if o.PerSlice
    fprintf('  series prefix    : "%s"  (must be WITHOUT slice timing correction)\n', prefix);
else
    fprintf('  series prefix    : "%s"  (slice-time corrected)\n', prefix);
end
fprintf('  TR               : %.4f s\n', TR);
fprintf('  design sampled at: bin %d of %d = %.4f s into each volume\n', T0, Tb, ref);
fprintf('  shift range      : %+.2f to %+.2f s, step %.5f s (%d points)\n', ...
        kRange(1)*dt, kRange(end)*dt, dt, numel(kRange));
fprintf('  events from      : %s\n', eventDir);

subs = list_dirs(derivRoot, 'sub-*', o.Subjects);
rows = {}; curves = {};

for si = 1:numel(subs)
    sub  = subs{si};
    sess = list_dirs(fullfile(derivRoot, sub), 'ses-*', o.Sessions);
    for se = 1:numel(sess)
        ses  = sess{se};
        fdir = fullfile(derivRoot, sub, ses, 'func');
        if ~isfolder(fdir), continue; end
        for tk = o.Tasks(:)'
            boldPath = find_one(fdir, sprintf('%ssub-*_%s_task-%d*_bold.nii', prefix, ses, tk));
            evPath   = find_one(eventDir, sprintf('%s_%s_task-%d*.txt', sub, ses, tk));
            rpPath   = find_one(fdir, sprintf('rp_*%s_task-%d*.txt', ses, tk));
            if isempty(boldPath) || isempty(evPath)
                fprintf('  skip %s %s task-%d  (bold found: %d, events found: %d)\n', ...
                        sub, ses, tk, ~isempty(boldPath), ~isempty(evPath));
                continue;
            end

            V   = spm_vol(boldPath);
            ns  = numel(V);
            dim = V(1).dim;
            Y   = reshape(single(spm_read_vols(V)), [], ns)';
            X   = design_matrix(evPath, dt, Tb, T0, ns, kRange);
            N   = nuisance(ns, TR, o.HPFSec, rpPath);

            if o.PerSlice
                st = slice_times(o.SliceTiming, dim(3));
                for z = 1:dim(3)
                    [dh, Fp, F0, nv] = sweep(Y(:, slice_idx(dim, z)), N, X, kRange, dt, o.KeepPct);
                    if nv == 0, continue; end
                    rows{end+1} = {sub, ses, tk, z, st(z), nv, dh, Fp, F0};   %#ok<AGROW>
                end
                fprintf('  %-11s %-7s task-%d  %d slices fitted\n', sub, ses, tk, dim(3));
            else
                [dh, Fp, F0, nv, curve] = sweep(Y, N, X, kRange, dt, o.KeepPct);
                if nv == 0
                    fprintf('  %-11s %-7s task-%d  no usable voxels\n', sub, ses, tk);
                else
                    rows{end+1}   = {sub, ses, tk, NaN, NaN, nv, dh, Fp, F0}; %#ok<AGROW>
                    curves{end+1} = curve;                                    %#ok<AGROW>
                    fprintf('  %-11s %-7s task-%d  n=%5d  delta = %+6.3f s  F peak %8.2f  (F at 0: %8.2f)\n', ...
                            sub, ses, tk, nv, dh, Fp, F0);
                end
            end
            clear Y X;
        end
    end
end

if isempty(rows)
    error('onset_check: nothing measured. Check derivRoot, eventDir and Prefix.');
end

T = cell2table(vertcat(rows{:}), 'VariableNames', ...
      {'subject_id','session','task','slice','slice_time_s','n_voxels', ...
       'delta_hat_s','F_peak','F_at_zero'});
writetable(T, csvPath);
fprintf('\n  wrote %s\n', csvPath);

if o.PerSlice
    ok = T.F_peak > 4;
    fprintf('\n  PER-SLICE FIT\n');
    if sum(ok) >= 6
        b   = [ones(sum(ok),1) T.slice_time_s(ok)] \ T.delta_hat_s(ok);
        tau = b(1) + ref;
        fprintf('    %d of %d slice-runs above F > 4\n', sum(ok), height(T));
        fprintf('    delta_hat = %+.3f  +  %.3f * SliceTiming\n', b(1), b(2));
        fprintf('    slope 1.00 would confirm SliceTiming describes the acquisition.\n');
        fprintf('    event-clock offset = intercept + %.3f = %+.3f s = %+.3f volumes\n', ...
                ref, tau, tau/TR);
    else
        fprintf('    skipped: only %d slice-runs reached F > 4.\n', sum(ok));
    end
else
    C = vertcat(curves{:});
    good = max(C, [], 2) > 0;
    C = C(good, :) ./ max(C(good, :), [], 2);
    m  = mean(C, 1);
    dh = refine(kRange, m, dt);
    try
        writematrix([kRange(:)*dt, m(:)], fullfile(outDir, 'onset_check_curve.csv'));
    catch
        csvwrite(fullfile(outDir, 'onset_check_curve.csv'), [kRange(:)*dt, m(:)]);
    end
    sd = NaN; if height(T) > 1, sd = std(T.delta_hat_s); end
    fprintf('\n  POOLED OVER %d RUNS\n', size(C,1));
    fprintf('    delta_hat = %+.3f s  =  %+.3f volumes\n', dh, dh/TR);
    fprintf('    per-run spread (sd) = %.3f s\n', sd);
    if abs(dh/TR) < 0.35
        fprintf('    VERDICT: event time 0 is the first saved volume. No correction needed.\n');
    else
        fprintf('    VERDICT: offset by %+.2f volumes.\n', dh/TR);
        fprintf('             Correct every onset as  onset - (%+.3f)  seconds.\n', dh);
    end
    fprintf('    curve written to %s\n', fullfile(outDir, 'onset_check_curve.csv'));
end
end

% ---------------------------------------------------------------------------

function X = design_matrix(evPath, dt, Tb, T0, ns, kRange)
% One HRF-convolved column per event name, pre-sampled at every shift.
% X is [ns x ncol x nshift]. The six-column read is GridCAT's own, so any
% event table GridCAT accepts is accepted here.
raw = fileread(evPath);
raw(raw == char(9)) = ',';
raw(raw == ';')     = ',';
C   = textscan(raw, '%s %f %f %f %f %f', 'Delimiter', ',');
% textscan can return a short last column, so trim all three to the same length.
nrow = min([numel(C{1}), numel(C{2}), numel(C{3})]);
nm   = C{1}(1:nrow);
on   = C{2}(1:nrow);
du   = C{3}(1:nrow);
keep = ~isnan(on) & ~isnan(du);
if ~any(keep), error('onset_check: no usable rows in %s', evPath); end
[names, ~, gid] = unique(nm(keep), 'stable');
on = on(keep); du = du(keep);

pad  = abs(kRange(1)) + Tb;
hrf  = spm_hrf(dt);
nbin = pad + (ns-1)*Tb + T0 + kRange(end) + numel(hrf) + Tb;
U    = zeros(nbin, numel(names));
for i = 1:numel(on)
    a = max(pad + round(on(i)/dt) + 1, 1);
    b = min(max(pad + round((on(i) + max(du(i), dt))/dt), a), nbin);
    U(a:b, gid(i)) = U(a:b, gid(i)) + 1;
end
for c = 1:size(U, 2)
    v = conv(U(:, c), hrf);
    U(:, c) = v(1:nbin);
end
X = zeros(ns, numel(names), numel(kRange));
for j = 1:numel(kRange)
    X(:, :, j) = U(pad + (0:ns-1)*Tb + T0 + kRange(j), :);
end
end

function N = nuisance(ns, TR, hpf, rpPath)
% Constant plus the DCT high-pass SPM would use, plus motion and its
% first differences. Identical for every shift, so it cannot bias the peak.
N = spm_dctmtx(ns, fix(2*(ns*TR)/hpf + 1));
if ~isempty(rpPath)
    rp = load(rpPath);
    rp = rp(1:min(ns, size(rp,1)), :);
    if size(rp,1) < ns, rp = [rp; repmat(rp(end,:), ns - size(rp,1), 1)]; end
    N = [N, rp, [zeros(1, size(rp,2)); diff(rp)]];
end
N = N(:, any(abs(N) > 0, 1));
end

function [dh, Fp, F0, nv, curve] = sweep(Y, N, X, kRange, dt, keepPct)
dh = NaN; Fp = NaN; F0 = NaN; nv = 0; curve = [];
mu  = mean(Y, 1);
pos = mu(mu > 0);
if isempty(pos), return; end
inb = mu > 0.2 * mean(pos);
if sum(inb) < 20, return; end

Y  = double(Y(:, inb));
PN = pinv(N);
R  = Y - N*(PN*Y);
sd = std(R, 0, 1);
ss = sort(sd);
thr = ss(max(1, min(numel(ss), round(numel(ss) * (1 - keepPct/100)))));
sel = sd >= thr;
if sum(sel) < 10, return; end
nv  = sum(sel);
yr  = mean(R(:, sel), 2);
syy = yr' * yr;
if syy <= 0, nv = 0; return; end

rN = rank(N);
curve = zeros(1, numel(kRange));
for j = 1:numel(kRange)
    Xr = X(:, :, j);
    Xr = Xr - N*(PN*Xr);
    q  = rank(Xr);
    if q == 0, continue; end
    r   = yr - Xr*(pinv(Xr)*yr);
    rss = r' * r;
    dfe = numel(yr) - q - rN;
    if rss <= 0 || dfe <= 0, continue; end
    curve(j) = ((syy - rss)/q) / (rss/dfe);
end
[Fp, i] = max(curve);
F0 = curve(kRange == 0);
dh = refine(kRange, curve, dt);
if i == 1 || i == numel(kRange)
    warning('onset_check:EdgeOfRange', ...
            'Peak sits at the edge of the shift range. Increase RangeTR.');
end
end

function dh = refine(kRange, curve, dt)
% Parabolic interpolation through the peak and its two neighbours.
[~, i] = max(curve);
if i > 1 && i < numel(curve)
    y1 = curve(i-1); y2 = curve(i); y3 = curve(i+1);
    den = y1 - 2*y2 + y3;
    s = 0;
    if den ~= 0, s = max(min(0.5*(y1 - y3)/den, 1), -1); end
    dh = (kRange(i) + s) * dt;
else
    dh = kRange(i) * dt;
end
end

function st = slice_times(given, nz)
if isempty(given)
    error(['onset_check: PerSlice needs SliceTiming - either a vector of %d ' ...
           'times in seconds or the path to a BOLD JSON sidecar.'], nz);
end
if ischar(given) || isstring(given)
    j = jsondecode(fileread(char(given)));
    if ~isfield(j, 'SliceTiming')
        error('onset_check: %s has no SliceTiming.', char(given));
    end
    given = j.SliceTiming;
end
st = given(:);
if numel(st) ~= nz
    error('onset_check: SliceTiming has %d entries but the image has %d slices.', ...
          numel(st), nz);
end
end

function idx = slice_idx(dim, z)
[a, b] = ndgrid(1:dim(1), 1:dim(2));
idx = sub2ind(dim, a(:), b(:), repmat(z, numel(a), 1));
end

function out = list_dirs(root, pat, keep)
d = dir(fullfile(root, pat));
if isempty(d), out = {}; return; end
out = sort({d([d.isdir]).name});
if ~isempty(keep), out = out(ismember(out, keep)); end
end

function f = find_one(dirPath, pat)
f = '';
if isempty(dirPath) || ~isfolder(dirPath), return; end
g = dir(fullfile(dirPath, pat));
if ~isempty(g), g = g(~startsWith({g.name}, '._')); end
if isempty(g), return; end
[~, k] = max([g.datenum]);
f = fullfile(g(k).folder, g(k).name);
end
