function write_run_manifest(DER, SUB, SES, CFG, boldPicks, others, anchor)
% WRITE_RUN_MANIFEST  Record which run of each scan this session used.
%
% Written as <ses>/preproc_run_selection.tsv. Two readers need it:
%   - prepare_gridcat_directory.m, to pair each task with its event table even
%     when the event file's name says nothing about runs;
%   - you, six months later, when a session looks odd and the question is
%     which T2w it was coregistered to.
%
% Rows are merged, not replaced. Running the coreg stage on its own decides
% the T2w and the ROI masks but nothing about the BOLD, and the record of
% which BOLD run the realignment used has to survive that.

outFile = fullfile(DER.sesDir, 'preproc_run_selection.tsv');

rows = read_existing_manifest(outFile);

for i = 1:numel(boldPicks)
    if isempty(boldPicks(i).path), continue; end
    key = ['task-' boldPicks(i).label];
    rows.(matlab.lang.makeValidName(key)) = sprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s', ...
        SUB, SES, key, boldPicks(i).runLabel, clock_of(boldPicks(i).acqSec), ...
        boldPicks(i).path, boldPicks(i).reason);
end

types = fieldnames(others);
for i = 1:numel(types)
    f = others.(types{i});
    if isempty(f), continue; end
    ent = bids_entities(f);
    runLabel = ent.run;
    if isempty(runLabel), runLabel = 'no-run'; end
    key = strrep(types{i}, '_', '-');
    rows.(matlab.lang.makeValidName(key)) = sprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s', ...
        SUB, SES, key, runLabel, clock_of(bids_acq_time(f)), f, ...
        'matched to the functional runs');
end

fid = fopen(outFile, 'w');
if fid < 0
    warning('run_spm_preproc:ManifestFailed', 'Could not write %s', outFile);
    return;
end
cleanupFid = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'subject\tsession\ttype\tselected_run\tacq_time\tsource_file\treason\n');
keys = sort(fieldnames(rows));
for i = 1:numel(keys)
    fprintf(fid, '%s\n', rows.(keys{i}));
end
fprintf(fid, '#\tanchor\t%s\t%s\tgap=%gmin\tstrict=%d\t\n', ...
    anchor.name, clock_of(anchor.sec), CFG.RUN_MATCH_GAP_SEC/60, CFG.STRICT_RUN_MATCHING);

fprintf('Run selection written to: %s\n', outFile);
end

function rows = read_existing_manifest(outFile)
% READ_EXISTING_MANIFEST  Previous rows, keyed by type, so they can be kept.
rows = struct();
if ~isfile(outFile), return; end

fid = fopen(outFile, 'r');
if fid == -1, return; end
cleanupFid = onCleanup(@() fclose(fid)); %#ok<NASGU>

lineNo = 0;
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    lineNo = lineNo + 1;
    if lineNo == 1, continue; end                    % header
    if isempty(strtrim(line)), continue; end
    if line(1) == '#', continue; end                 % the anchor footer

    parts = strsplit(line, sprintf('\t'));
    if numel(parts) < 3, continue; end
    key = matlab.lang.makeValidName(strtrim(parts{3}));
    rows.(key) = line;
end
end
