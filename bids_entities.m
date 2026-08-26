function E = bids_entities(fileName)
% BIDS_ENTITIES  Split a BIDS filename into its entities.
%
%   E = bids_entities('sub-01s13_ses-01_task-run1_run-2_bold.nii.gz')
%
%   E.name     'sub-01s13_ses-01_task-run1_run-2_bold.nii.gz'
%   E.stem     'sub-01s13_ses-01_task-run1_run-2_bold'
%   E.ext      '.nii.gz'
%   E.suffix   'bold'
%   E.sub      'sub-01s13'      (with the prefix, as it appears in paths)
%   E.ses      'ses-01'
%   E.task     'run1'           (label only, without 'task-')
%   E.run      'run-2'          ('' when the file carries no run entity)
%   E.runNum   2                (NaN when there is no run entity)
%   E.acq      ''               likewise for acq-, dir-, echo-, part-
%   E.prefix   ''               anything prepended by SPM, e.g. 'su' or 'rp_'
%
%   Anything before 'sub-' is treated as a prefix, so a file SPM renamed
%   ('su' + the BIDS name, 'rp_' + the BIDS name) still parses and E.prefix
%   tells the caller it is derived rather than raw.
%
%   See also BIDS_LIST_RUNS, BIDS_PICK_RUN, BIDS_NAME.

name = basename(char(fileName));

% Strip the extension, keeping .nii.gz together
ext = '';
stem = name;
lowerName = lower(name);
knownExt = {'.nii.gz', '.nii', '.json', '.txt', '.tsv', '.img', '.hdr', '.mat'};
for i = 1:numel(knownExt)
    if endsWith(lowerName, knownExt{i})
        ext = name(end-numel(knownExt{i})+1:end);
        stem = name(1:end-numel(knownExt{i}));
        break;
    end
end

E = struct('name', name, 'stem', stem, 'ext', ext, 'suffix', '', 'prefix', '', ...
           'sub', '', 'ses', '', 'task', '', 'run', '', 'runNum', NaN, ...
           'acq', '', 'dir', '', 'echo', '', 'part', '');

% Anything in front of 'sub-' is a prefix SPM (or this pipeline) added
bidsPart = stem;
subStart = regexp(stem, 'sub-', 'once');
if ~isempty(subStart)
    if subStart > 1
        E.prefix = stem(1:subStart-1);
    end
    bidsPart = stem(subStart:end);
end

parts = strsplit(bidsPart, '_');

for i = 1:numel(parts)
    p = parts{i};
    if isempty(p), continue; end

    tok = regexp(p, '^([A-Za-z0-9]+)-(.+)$', 'tokens', 'once');
    if isempty(tok)
        % No key-value pair: the last such chunk is the BIDS suffix
        E.suffix = p;
        continue;
    end

    key = lower(tok{1});
    val = tok{2};
    switch key
        case 'sub',  E.sub  = ['sub-' val];
        case 'ses',  E.ses  = ['ses-' val];
        case 'task', E.task = val;
        case 'run'
            E.run = ['run-' val];
            n = str2double(val);
            if ~isnan(n), E.runNum = n; end
        case 'acq',  E.acq  = val;
        case 'dir',  E.dir  = val;
        case 'echo', E.echo = val;
        case 'part', E.part = val;
        otherwise
            % Unknown entities are ignored on purpose: the pipeline only ever
            % matches on the ones above.
    end
end
end

function b = basename(p)
% Basename that keeps the full extension (fileparts drops only the last one).
idx = find(p == filesep | p == '/', 1, 'last');
if isempty(idx)
    b = p;
else
    b = p(idx+1:end);
end
end
