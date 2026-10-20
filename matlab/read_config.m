function C = read_config(cfgFile)
% READ_CONFIG  Read config.sh into a struct of strings.
%
%   C = read_config('/path/to/config.sh')
%
% Every value comes back as a char array, exactly as written in the file.
% Nothing is renamed and nothing is guessed, so you always know what you
% have. Convert at the point of use:
%
%   str2double(C.RUN_MATCH_GAP_MIN)      a number
%   strcmp(C.STRICT_RUN_MATCHING,'true') a flag
%   strsplit(C.TASKS, ',')               a list
%
% SPM and GridCAT parameters are not here. They live in spm_options.m as
% real MATLAB values, so they need no parsing and no conversion.

if nargin < 1 || isempty(cfgFile)
    cfgFile = fullfile(fileparts(mfilename('fullpath')), 'config.sh');
end
assert(isfile(cfgFile), 'read_config: not found: %s', cfgFile);

C   = struct();
fid = fopen(cfgFile, 'r');
assert(fid ~= -1, 'read_config: cannot open: %s', cfgFile);
closeFile = onCleanup(@() fclose(fid));

while true
    line = fgetl(fid);
    if ~ischar(line), break; end

    line = strtrim(line);
    if isempty(line) || line(1) == '#'
        continue;
    end

    eq = find(line == '=', 1, 'first');     % first '=' only, so paths survive
    if isempty(eq), continue; end

    key = strtrim(line(1:eq-1));
    if ~isvarname(key), continue; end

    C.(key) = strtrim(line(eq+1:end));
end

end
