function js = read_json_file(jsonFile)
% READ_JSON_FILE  Read a JSON sidecar into a struct.
%
%   Returns an empty struct when the file is missing or cannot be parsed, so
%   callers can check with isfield() instead of wrapping every read in a try.

js = struct();

if isempty(jsonFile) || ~isfile(jsonFile)
    return;
end

fid = fopen(jsonFile, 'r');
if fid == -1
    return;
end
raw = fread(fid, inf, '*char')';
fclose(fid);

try
    js = jsondecode(raw);
catch
    warning('read_json_file:ParseFailed', 'Could not parse JSON: %s', jsonFile);
    js = struct();
end
end
