function didReset = bids_reset_deriv(D, doReset)
% BIDS_RESET_DERIV  Make this session's preprocessing derivatives ready to use.
%
%   didReset = bids_reset_deriv(D, doReset)
%
%   D        from bids_deriv
%   doReset  true  -> delete everything under <preproc>/sub-XX/ses-YY first
%            false -> keep what is there (a resumed run reads it back)
%
%   A full preprocessing run always starts from an empty session directory, so
%   the output can never be a mixture of this run and the last one: no stale
%   voxel displacement maps, no u*_bold.nii from a different chosen run, no
%   left-over topup fields. Partial runs (--stages realign,coreg,smooth and
%   friends) keep the directory, because that is exactly the output they
%   resume from.
%
%   Only the session's own directory is removed, never the whole dataset, so
%   the array tasks of one submission cannot delete each other's work.
%
%   See also BIDS_DERIV, BIDS_STAGE.

didReset = false;

assert_safe_target(D);

if doReset && isfolder(D.sesDir)
    fprintf('Clearing previous derivatives: %s\n', D.sesDir);
    [ok, msg] = rmdir(D.sesDir, 's');
    if ~ok
        error('bids_reset_deriv:CannotRemove', ...
            'Could not clear %s\n  %s', D.sesDir, msg);
    end
    didReset = true;
end

for d = {D.sesDir, D.func, D.anat, D.fmap}
    if ~isfolder(d{1})
        [ok, msg] = mkdir(d{1});
        if ~ok
            error('bids_reset_deriv:CannotCreate', 'Could not create %s\n  %s', d{1}, msg);
        end
    end
end

write_dataset_description(D);
end

% ============================== helpers ==============================

function assert_safe_target(D)
% A mistyped DERIV_ROOT must not turn into "rm -rf" on something valuable.
sesDir = strip_trailing_sep(D.sesDir);
preproc = strip_trailing_sep(D.preproc);

if isempty(D.sub)
    error('bids_reset_deriv:NoSubject', ...
        'Refusing to touch %s without a subject — bids_deriv was called without SUB/SES.', preproc);
end
if ~startsWith(sesDir, preproc)
    error('bids_reset_deriv:OutsideRoot', ...
        'Refusing to clear %s: it is not inside the derivatives dataset %s', sesDir, preproc);
end
% Guard against DERIV_ROOT=/ or similar
if numel(strsplit(strrep(sesDir, '\', '/'), '/')) < 4
    error('bids_reset_deriv:TooShallow', ...
        'Refusing to clear a top-level path: %s\nCheck DERIV_ROOT in pipeline_config.cfg.', sesDir);
end
end

function s = strip_trailing_sep(s)
while ~isempty(s) && (s(end) == '/' || s(end) == filesep)
    s(end) = [];
end
end

function write_dataset_description(D)
% A BIDS derivatives dataset is supposed to say what made it. Written once;
% cheap enough to rewrite whenever a session is prepared.
f = fullfile(D.preproc, 'dataset_description.json');
if isfile(f), return; end

txt = sprintf([ ...
'{\n' ...
'  "Name": "%s",\n' ...
'  "BIDSVersion": "1.8.0",\n' ...
'  "DatasetType": "derivative",\n' ...
'  "Description": "SPM preprocessing: distortion correction, realignment, coregistration, smoothing.",\n' ...
'  "GeneratedBy": [\n' ...
'    {\n' ...
'      "Name": "Grid-Cell-Pipeline",\n' ...
'      "Description": "run_spm_preproc.m",\n' ...
'      "CodeURL": "https://github.com/ju-li-lab/Grid-Cell-Pipeline"\n' ...
'    }\n' ...
'  ]\n' ...
'}\n'], last_folder(D.preproc));

fid = fopen(f, 'w');
if fid > 0
    fprintf(fid, '%s', txt);
    fclose(fid);
end
end

function n = last_folder(p)
p = strip_trailing_sep(p);
parts = strsplit(strrep(p, '\', '/'), '/');
n = parts{end};
end
