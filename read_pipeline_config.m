function cfg = read_pipeline_config(cfgFile)
% READ_PIPELINE_CONFIG  Read the central pipeline configuration file.
%
% Usage:
%   cfg = read_pipeline_config('pipeline_config.cfg')
%   cfg = read_pipeline_config('/full/path/to/pipeline_config.cfg')
%
% Returns:
%   cfg - A struct where each KEY becomes a field name.
%         Numeric values are automatically converted to numbers (int or float).
%         Comma-separated values become numeric arrays: "1,2,3" -> [1 2 3]
%         String values remain as char arrays.
%         Empty values are stored as empty strings: ''
%
% Examples:
%   cfg = read_pipeline_config('pipeline_config.cfg');
%   cfg.BIDS_ROOT         % '/sc-projects/...' (string)
%   cfg.TASKS             % [1 2 3] (numeric array)
%   cfg.TR                % 2 (numeric scalar)
%   cfg.SPM_DIR           % '/sc-projects/...' (string)
%   cfg.MASKING_THRESHOLD % 0.8 (numeric scalar, float)
%
% Features:
%   - Ignores comment lines (starting with #)
%   - Handles empty lines and blank lines gracefully
%   - Strips whitespace from keys and values
%   - Auto-detects and converts numeric values (integers and floats)
%   - Converts comma-separated numbers to numeric arrays
%   - Converts comma-separated non-numeric values to cell arrays of strings
%   - Keeps string values as char arrays
%   - Handles paths with '=' signs (only splits on first '=')
%   - Provides clear error messages if file is missing or malformed
%
% Note:
%   Logical values like 'true' or 'false' remain as strings.
%   To use them in MATLAB, check: strcmp(cfg.DRY_RUN, 'true')

    % Check if file was provided
    if nargin < 1 || isempty(cfgFile)
        error('read_pipeline_config:NoFile', ...
              'Configuration file path required.\nUsage: cfg = read_pipeline_config(cfgFile)');
    end

    % Check if file exists
    if ~isfile(cfgFile)
        error('read_pipeline_config:FileNotFound', ...
              sprintf('Configuration file not found:\n  %s\n\nCheck the path and try again.', cfgFile));
    end

    % Initialize output struct
    cfg = struct();

    try
        % Read file line by line
        fid = fopen(cfgFile, 'r');
        if fid == -1
            error('read_pipeline_config:CannotOpen', ...
                  sprintf('Cannot open configuration file:\n  %s', cfgFile));
        end

        lineNum = 0;

        while ~feof(fid)
            lineNum = lineNum + 1;
            line = fgetl(fid);

            % Handle end of file
            if ~ischar(line)
                break;
            end

            % Strip leading/trailing whitespace
            line = strtrim(line);

            % Skip empty lines and comments
            if isempty(line) || line(1) == '#'
                continue;
            end

            % Find the first '=' to split key and value
            eqPos = find(line == '=', 1, 'first');

            if isempty(eqPos)
                % Line has no '=' sign, skip with warning
                warning('read_pipeline_config:MalformedLine', ...
                        sprintf('Line %d has no ''='' sign, skipping:\n  %s', lineNum, line));
                continue;
            end

            % Extract key and value
            key = strtrim(line(1:eqPos-1));
            val = strtrim(line(eqPos+1:end));

            % Validate key (must be non-empty and valid MATLAB identifier)
            if isempty(key)
                warning('read_pipeline_config:EmptyKey', ...
                        sprintf('Line %d has empty key, skipping.', lineNum));
                continue;
            end

            % Check if key is a valid MATLAB identifier
            if ~isvarname(key)
                warning('read_pipeline_config:InvalidKey', ...
                        sprintf('Line %d: ''%s'' is not a valid MATLAB identifier, skipping.', ...
                                lineNum, key));
                continue;
            end

            % Try to convert value to numeric type
            convertedVal = tryConvertToNumeric(val);

            % Store in struct
            cfg.(key) = convertedVal;
        end

        fclose(fid);

    catch ME
        if ~strcmp(ME.identifier, 'read_pipeline_config:FileNotFound') && ...
           ~strcmp(ME.identifier, 'read_pipeline_config:CannotOpen') && ...
           ~strcmp(ME.identifier, 'read_pipeline_config:NoFile')
            fclose(fid);
        end
        rethrow(ME);
    end

    if isempty(fieldnames(cfg))
        warning('read_pipeline_config:EmptyConfig', ...
                'Configuration file appears to be empty or contains only comments.');
    end

end

% ========================================================================
% Helper function: Try to convert a string value to numeric type
% ========================================================================

function val = tryConvertToNumeric(strVal)
% TRYCONVERTONUMERIC  Attempt numeric conversion of string value.
%
% Handles:
%   - Single numbers: "2" → 2
%   - Floats: "0.8" → 0.8
%   - Comma-separated numbers: "1,2,3" → [1 2 3]
%   - Returns original string if conversion fails

    if ~ischar(strVal) && ~isstring(strVal)
        val = strVal;
        return;
    end

    % Convert to char if string
    strVal = char(strVal);

    % Try comma-separated numeric array first
    if contains(strVal, ',')
        parts = strsplit(strtrim(strVal), ',');
        numArray = [];
        allNumeric = true;

        for i = 1:length(parts)
            part = strtrim(parts{i});
            num = str2double(part);
            if isnan(num)
                allNumeric = false;
                break;
            end
            numArray = [numArray, num]; %#ok<AGROW>
        end

        if allNumeric && ~isempty(numArray)
            val = numArray;
            return;
        end
    end

    % Try single numeric value
    num = str2double(strVal);
    if ~isnan(num)
        val = num;
    else
        % Keep as string
        val = strVal;
    end

end
