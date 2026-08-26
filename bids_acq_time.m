function [acqSec, series] = bids_acq_time(niiFile)
% BIDS_ACQ_TIME  When a scan was acquired, from its JSON sidecar.
%
%   [acqSec, series] = bids_acq_time('/path/sub-01_ses-01_run-2_T2w.nii')
%
%   acqSec  seconds since midnight (AcquisitionDateTime or AcquisitionTime)
%   series  the scanner's SeriesNumber
%
%   Both are NaN when the sidecar is missing or does not carry the field.
%
%   This is what lets the pipeline tell apart runs acquired back to back from
%   runs acquired after the subject climbed out of the scanner and came back:
%   the second block starts minutes to hours later. Run numbers alone cannot
%   express that — run-2 of the T2w and run-2 of a task may well belong to
%   different blocks — so acquisition time is the thing that matches them up.
%
%   See also BIDS_LIST_RUNS, BIDS_PICK_RUN.

acqSec = NaN;
series = NaN;

if isempty(niiFile), return; end

jsonFile = regexprep(char(niiFile), '\.nii(\.gz)?$', '.json');
if ~isfile(jsonFile), return; end

js = read_json_file(jsonFile);
if isempty(fieldnames(js)), return; end

% Time of day: 'HH:MM:SS.ffffff', or an ISO datetime ending in one
raw = '';
if isfield(js, 'AcquisitionDateTime') && ~isempty(js.AcquisitionDateTime)
    raw = js.AcquisitionDateTime;
elseif isfield(js, 'AcquisitionTime') && ~isempty(js.AcquisitionTime)
    raw = js.AcquisitionTime;
end

if ~isempty(raw)
    raw = strtrim(char(raw));
    tok = regexp(raw, '(\d{2}):(\d{2}):(\d{2}(?:\.\d+)?)', 'tokens', 'once');
    if isempty(tok)
        % Bare DICOM form: HHMMSS.ffffff
        tok = regexp(raw, '^(\d{2})(\d{2})(\d{2}(?:\.\d+)?)$', 'tokens', 'once');
    end
    if ~isempty(tok)
        acqSec = str2double(tok{1}) * 3600 + str2double(tok{2}) * 60 + str2double(tok{3});
    end
end

if isfield(js, 'SeriesNumber') && ~isempty(js.SeriesNumber)
    v = js.SeriesNumber;
    if ischar(v), v = str2double(v); end
    if isnumeric(v) && isscalar(v), series = double(v); end
end
end
