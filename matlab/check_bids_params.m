function check_bids_params(boldFiles, reverseFile, phasediffFile, fieldmapFile, CFG)
% CHECK_BIDS_PARAMS  Compare the config against the sidecars of the files
% this run actually chose.
%
% The files are the staged copies in the derivatives, and their JSON sidecars
% came along with them. Checking those rather than "whichever sidecar turns up
% first in the folder" matters as soon as a session has several runs: the run
% being preprocessed is the one whose acquisition parameters have to match the
% config, and a mismatch used to hide behind another run's sidecar.
%
% Nothing here stops the pipeline. The config values are always the ones used,
% because they are sometimes deliberate overrides.

% Acquisition parameters live in spm_options.m, in the shape SPM wants them.
fm = CFG.spm.fieldmap.defaultsval;   % .et = [TE1 TE2], .tert, .blipdir

fprintf('\n--- Validating config against the BIDS JSON sidecars ---\n');
nWarnings = 0;

% --- The BOLD run being preprocessed ---
boldJson = '';
for i = 1:numel(boldFiles)
    if isempty(boldFiles{i}), continue; end
    cand = regexprep(boldFiles{i}, '\.nii(\.gz)?$', '.json');
    if isfile(cand)
        boldJson = cand;
        break;
    end
end

bidsPE = '';
bidsPE_hasSign = false;

if ~isempty(boldJson)
    bj = read_json_file(boldJson);
    fprintf('  BOLD JSON: %s\n', boldJson);

    % -- TR --
    if isfield(bj, 'RepetitionTime')
        bidsVal = bj.RepetitionTime;
        if isfield(CFG, 'TR') && abs(CFG.TR - bidsVal) > 0.001
            fprintf('  WARNING: TR mismatch — config: %.4f s, BIDS JSON: %.4f s\n', CFG.TR, bidsVal);
            nWarnings = nWarnings + 1;
        else
            fprintf('  TR: %.4f s (matches)\n', bidsVal);
        end
    end

    % -- Phase encoding direction --
    % BIDS has two PE fields:
    %   PhaseEncodingDirection = full direction with sign (i, i-, j, j-, k, k-)
    %   PhaseEncodingAxis      = axis only, no sign (i, j, k)
    % We use Direction if available, otherwise fall back to Axis.
    if isfield(bj, 'PhaseEncodingDirection') && ~isempty(bj.PhaseEncodingDirection)
        bidsPE = bj.PhaseEncodingDirection;
        bidsPE_hasSign = true;
    elseif isfield(bj, 'PhaseEncodingAxis') && ~isempty(bj.PhaseEncodingAxis)
        bidsPE = bj.PhaseEncodingAxis;
        bidsPE_hasSign = false;
    end
    if ~isempty(bidsPE)
        fslPE = bids_pe_to_fsl(bidsPE);
        bidsAxis = regexprep(fslPE, '-$', '');

        if strcmp(CFG.PREPROC_MODE, 'topup') && ~isempty(fslPE)
            configAxis = regexprep(CFG.TOPUP_PE_DIR_BOLD, '-$', '');  % 'y-' -> 'y'
            if bidsPE_hasSign
                if ~strcmp(fslPE, CFG.TOPUP_PE_DIR_BOLD)
                    fprintf('  WARNING: PE direction mismatch — config TOPUP_PE_DIR_BOLD: %s, BIDS JSON: %s (=%s in FSL)\n', ...
                        CFG.TOPUP_PE_DIR_BOLD, bidsPE, fslPE);
                    nWarnings = nWarnings + 1;
                else
                    fprintf('  PE direction (BOLD): %s = %s in FSL (matches)\n', bidsPE, fslPE);
                end
            else
                if ~strcmp(bidsAxis, configAxis)
                    fprintf('  WARNING: PE axis mismatch — config TOPUP_PE_DIR_BOLD axis: %s, BIDS JSON: %s (=%s in FSL)\n', ...
                        configAxis, bidsPE, bidsAxis);
                    nWarnings = nWarnings + 1;
                else
                    fprintf('  PE axis (BOLD): %s (matches). Note: BIDS JSON has PhaseEncodingAxis only (no sign) — cannot verify +/- direction.\n', bidsAxis);
                end
            end
        end

        if strcmp(CFG.PREPROC_MODE, 'realign_unwarp') && ~isempty(fslPE)
            if bidsPE_hasSign
                % Check blip direction consistency: j -> +1, j- -> -1
                expectedBlip = 1;
                if endsWith(bidsPE, '-'), expectedBlip = -1; end
                if fm.blipdir ~= expectedBlip
                    fprintf('  WARNING: BLIP_DIRECTION mismatch — config: %d, BIDS JSON PE=%s suggests: %d\n', ...
                        fm.blipdir, bidsPE, expectedBlip);
                    nWarnings = nWarnings + 1;
                else
                    fprintf('  BLIP_DIRECTION: %d (consistent with BIDS PE=%s)\n', fm.blipdir, bidsPE);
                end
            else
                fprintf('  BIDS JSON has PhaseEncodingAxis=%s only (no sign) — cannot verify BLIP_DIRECTION.\n', bidsPE);
            end
        end
    end

    % -- Total readout time (topup) --
    if strcmp(CFG.PREPROC_MODE, 'topup')
        bidsReadout = [];
        if isfield(bj, 'EstimatedTotalReadoutTime')
            bidsReadout = bj.EstimatedTotalReadoutTime;
        elseif isfield(bj, 'TotalReadoutTime')
            bidsReadout = bj.TotalReadoutTime;
        end
        if ~isempty(bidsReadout)
            if abs(CFG.TOPUP_READOUT_SEC - bidsReadout) > 0.001
                fprintf('  WARNING: Readout time mismatch — config TOPUP_READOUT_SEC: %.6f s, BIDS JSON: %.6f s\n', ...
                    CFG.TOPUP_READOUT_SEC, bidsReadout);
                nWarnings = nWarnings + 1;
            else
                fprintf('  Readout time (topup): %.6f s (matches)\n', bidsReadout);
            end
            % Sanity: warn if value looks like milliseconds
            if CFG.TOPUP_READOUT_SEC > 1
                fprintf('  WARNING: TOPUP_READOUT_SEC=%.3f looks like milliseconds — should be in seconds\n', ...
                    CFG.TOPUP_READOUT_SEC);
                nWarnings = nWarnings + 1;
            end
        end
    end

    % -- Total readout time (realign_unwarp — in ms) --
    if strcmp(CFG.PREPROC_MODE, 'realign_unwarp')
        bidsReadout = [];
        if isfield(bj, 'TotalReadoutTime')
            bidsReadout = bj.TotalReadoutTime * 1000;  % BIDS stores seconds
        elseif isfield(bj, 'EstimatedTotalReadoutTime')
            bidsReadout = bj.EstimatedTotalReadoutTime * 1000;
        end
        if ~isempty(bidsReadout)
            if abs(fm.tert - bidsReadout) > 1  % tolerance 1 ms
                fprintf('  WARNING: Readout time mismatch — config TOTAL_READOUT_MS: %.3f ms, BIDS JSON: %.3f ms\n', ...
                    fm.tert, bidsReadout);
                nWarnings = nWarnings + 1;
            else
                fprintf('  Readout time (VDM): %.3f ms (matches)\n', bidsReadout);
            end
        end
    end
else
    fprintf('  No BOLD JSON sidecar found — skipping BOLD parameter validation.\n');
end

% --- The reverse-PE EPI topup was given ---
if strcmp(CFG.PREPROC_MODE, 'topup') && ~isempty(reverseFile)
    revJson = regexprep(reverseFile, '\.nii(\.gz)?$', '.json');

    if isfile(revJson)
        rj = read_json_file(revJson);
        fprintf('  Reverse-PE JSON: %s\n', revJson);

        revPE = '';
        revPE_hasSign = false;
        if isfield(rj, 'PhaseEncodingDirection') && ~isempty(rj.PhaseEncodingDirection)
            revPE = rj.PhaseEncodingDirection;
            revPE_hasSign = true;
        elseif isfield(rj, 'PhaseEncodingAxis') && ~isempty(rj.PhaseEncodingAxis)
            revPE = rj.PhaseEncodingAxis;
            revPE_hasSign = false;
        end
        if ~isempty(revPE)
            fslRevPE = bids_pe_to_fsl(revPE);
            configRevAxis = regexprep(CFG.TOPUP_PE_DIR_REVERSE, '-$', '');
            bidsRevAxis = regexprep(fslRevPE, '-$', '');

            if revPE_hasSign
                if ~isempty(fslRevPE) && ~strcmp(fslRevPE, CFG.TOPUP_PE_DIR_REVERSE)
                    fprintf('  WARNING: Reverse PE direction mismatch — config: %s, BIDS JSON: %s (=%s in FSL)\n', ...
                        CFG.TOPUP_PE_DIR_REVERSE, revPE, fslRevPE);
                    nWarnings = nWarnings + 1;
                else
                    fprintf('  PE direction (reverse): %s = %s in FSL (matches)\n', revPE, fslRevPE);
                end
            else
                if ~strcmp(bidsRevAxis, configRevAxis)
                    fprintf('  WARNING: Reverse PE axis mismatch — config axis: %s, BIDS JSON: %s (=%s in FSL)\n', ...
                        configRevAxis, revPE, bidsRevAxis);
                    nWarnings = nWarnings + 1;
                else
                    fprintf('  PE axis (reverse): %s (matches). PhaseEncodingAxis only — cannot verify +/- sign.\n', bidsRevAxis);
                end
            end
        end

        % Forward and reverse must be opposite, or topup has nothing to work with
        if bidsPE_hasSign && revPE_hasSign && ~isempty(bidsPE) && ~isempty(revPE)
            fwdVec = pe_dir_to_vector(bids_pe_to_fsl(bidsPE));
            revVec = pe_dir_to_vector(bids_pe_to_fsl(revPE));
            if ~all(fwdVec + revVec == 0)
                fprintf('  WARNING: Forward PE (%s) and reverse PE (%s) are NOT opposite!\n', bidsPE, revPE);
                fprintf('           Forward vector: [%d %d %d], Reverse vector: [%d %d %d]\n', ...
                    fwdVec, revVec);
                nWarnings = nWarnings + 1;
            else
                fprintf('  Forward/reverse PE vectors are opposite (OK)\n');
            end
        elseif ~isempty(bidsPE) && ~isempty(revPE)
            fwdAxis = regexprep(bids_pe_to_fsl(bidsPE), '-$', '');
            revAxis = regexprep(bids_pe_to_fsl(revPE), '-$', '');
            if strcmp(fwdAxis, revAxis)
                fprintf('  Forward and reverse PE share axis %s (OK). Only PhaseEncodingAxis available — verify opposite signs manually.\n', fwdAxis);
            else
                fprintf('  WARNING: Forward PE axis (%s) and reverse PE axis (%s) are on DIFFERENT axes!\n', fwdAxis, revAxis);
                nWarnings = nWarnings + 1;
            end
        end

        % The two EPIs must also come from the same visit to the scanner, or
        % topup models the head moving as a distortion field.
        nWarnings = nWarnings + warn_if_unmatched(boldFiles, reverseFile, 'reverse-PE EPI', CFG);
    else
        fprintf('  No reverse-PE JSON sidecar found — skipping reverse-PE validation.\n');
    end
end

% --- The GRE phasediff ---
if strcmp(CFG.PREPROC_MODE, 'realign_unwarp') && ~isempty(phasediffFile)
    phasediffJson = regexprep(phasediffFile, '\.nii(\.gz)?$', '.json');
    if isfile(phasediffJson)
        pj = read_json_file(phasediffJson);
        fprintf('  Phasediff JSON: %s\n', phasediffJson);

        if isfield(pj, 'EchoTime1') && isfield(pj, 'EchoTime2')
            bidsTE1 = pj.EchoTime1 * 1000;  % to ms
            bidsTE2 = pj.EchoTime2 * 1000;
            if abs(fm.et(1) - bidsTE1) > 0.01
                fprintf('  WARNING: TE1 mismatch — config TE_SHORT_MS: %.3f ms, BIDS JSON EchoTime1: %.3f ms\n', ...
                    fm.et(1), bidsTE1);
                nWarnings = nWarnings + 1;
            else
                fprintf('  TE1: %.3f ms (matches)\n', bidsTE1);
            end
            if abs(fm.et(2) - bidsTE2) > 0.01
                fprintf('  WARNING: TE2 mismatch — config TE_LONG_MS: %.3f ms, BIDS JSON EchoTime2: %.3f ms\n', ...
                    fm.et(2), bidsTE2);
                nWarnings = nWarnings + 1;
            else
                fprintf('  TE2: %.3f ms (matches)\n', bidsTE2);
            end
        else
            fprintf('  Phasediff JSON has no EchoTime1/EchoTime2 — skipping TE validation.\n');
        end
    else
        fprintf('  No phasediff JSON sidecar — skipping TE validation.\n');
    end
    nWarnings = nWarnings + warn_if_unmatched(boldFiles, phasediffFile, 'phasediff', CFG);
end

% --- The precalculated fieldmap ---
if strcmp(CFG.PREPROC_MODE, 'precalc_fieldmap') && ~isempty(fieldmapFile)
    fmJson = regexprep(fieldmapFile, '\.nii(\.gz)?$', '.json');
    if ~isfile(fmJson)
        fprintf('  No fieldmap JSON next to %s — cannot cross-check units.\n', fieldmapFile);
        fprintf('  Using FIELDMAP_UNITS from config: %s\n', CFG.FIELDMAP_UNITS);
    else
        fj = read_json_file(fmJson);
        fprintf('  Fieldmap JSON: %s\n', fmJson);
        if isfield(fj, 'Units') && ~isempty(fj.Units)
            jsonUnits = lower(strtrim(fj.Units));
            cfgUnits  = CFG.FIELDMAP_UNITS;
            % 'rad/s' and 'rads' mean the same thing here
            normalise = @(u) regexprep(lower(u), '^(rad/s|rads|rad_per_s)$', 'rad/s');
            if ~strcmp(normalise(jsonUnits), normalise(cfgUnits))
                fprintf(['  WARNING: fieldmap units mismatch — config FIELDMAP_UNITS: %s, ' ...
                         'JSON Units: %s\n'], cfgUnits, jsonUnits);
                fprintf('           A wrong unit scales the whole distortion correction by 2*pi.\n');
                nWarnings = nWarnings + 1;
            else
                fprintf('  Fieldmap units: %s (matches config)\n', jsonUnits);
            end
        else
            fprintf('  Fieldmap JSON has no "Units" field — using config: %s\n', CFG.FIELDMAP_UNITS);
        end
        if isfield(fj, 'GeneratedBy') && isstruct(fj.GeneratedBy) && ...
           isfield(fj.GeneratedBy, 'MagnitudeSource')
            fprintf('  Magnitude source: %s\n', fj.GeneratedBy.MagnitudeSource);
        end
    end

    % A GRE-derived fieldmap is not EPI-based; getting this backwards inverts
    % how SPM matches the VDM to the EPI.
    if fm.epifm ~= 0
        fprintf(['  WARNING: EPI_BASED_FIELDMAP=%d with PREPROC_MODE=precalc_fieldmap.\n' ...
                 '           A fieldmap made from a GRE phasediff (make_fieldmaps.sh)\n' ...
                 '           is NOT EPI-based — set EPI_BASED_FIELDMAP=0.\n'], ...
                fm.epifm);
        nWarnings = nWarnings + 1;
    end
end

% --- Summary ---
if nWarnings == 0
    fprintf('  All BIDS JSON checks passed.\n');
else
    fprintf('\n  *** %d BIDS JSON warning(s) detected. ***\n', nWarnings);
    fprintf('  Config values will still be used. Check warnings above.\n');
    fprintf('  If the config values are intentional overrides, you can ignore these warnings.\n');
end
fprintf('\n');
end

function n = warn_if_unmatched(boldFiles, otherFile, label, CFG)
% WARN_IF_UNMATCHED  Say so when a chosen input was acquired far from the BOLD.
%
% The run selection already tries to match by acquisition time, but an explicit
% choice in run_selection.tsv overrides that, and a missing sidecar hides it.
% This is the last place to notice that a distortion correction is about to be
% estimated between scans from two different head positions.

n = 0;
if isempty(otherFile), return; end

otherSec = bids_acq_time(otherFile);
if isnan(otherSec), return; end

boldSecs = [];
for i = 1:numel(boldFiles)
    if isempty(boldFiles{i}), continue; end
    s = bids_acq_time(boldFiles{i});
    if ~isnan(s), boldSecs(end+1) = s; end %#ok<AGROW>
end
if isempty(boldSecs), return; end

gap = min(abs(boldSecs - otherSec));
if gap > CFG.RUN_MATCH_GAP_SEC
    fprintf(['  WARNING: the %s was acquired %.0f min from the nearest BOLD run.\n' ...
             '           They are probably from different visits to the scanner —\n' ...
             '           check run_selection.tsv.\n'], label, gap/60);
    n = 1;
end
end

function fslDir = bids_pe_to_fsl(bidsPE)
% BIDS_PE_TO_FSL  Convert BIDS PhaseEncodingDirection to FSL convention.
%   'i' -> 'x',  'i-' -> 'x-'
%   'j' -> 'y',  'j-' -> 'y-'
%   'k' -> 'z',  'k-' -> 'z-'
%   Already FSL format (x/y/z) passes through unchanged.
%   Returns '' if unknown.

switch bidsPE
    case 'i',  fslDir = 'x';
    case 'i-', fslDir = 'x-';
    case 'j',  fslDir = 'y';
    case 'j-', fslDir = 'y-';
    case 'k',  fslDir = 'z';
    case 'k-', fslDir = 'z-';
    case {'x','x-','y','y-','z','z-'}
        fslDir = bidsPE;  % Already FSL format
    otherwise
        fslDir = '';
        fprintf('  (Unknown BIDS PE direction: %s)\n', bidsPE);
end
end

% ===================== PROVENANCE LOG =====================
