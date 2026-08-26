function D = bids_deriv(cfgRaw, SUB, SES)
% BIDS_DERIV  Where this session's derived files live.
%
%   D = bids_deriv(cfgRaw, 'sub-01s13', 'ses-01')
%
%   cfgRaw is the struct from read_pipeline_config. Returns:
%
%     D.root        the derivatives root      <OUTPUT_ROOT>/derivatives
%     D.preproc     the preprocessing dataset <root>/spm-preproc
%     D.rois        imported ROI masks        <root>/rois
%     D.fieldmaps   fieldmaps built from structurals   <root>/fieldmaps
%
%     D.sesDir  <preproc>/sub-XX/ses-YY
%     D.func    <preproc>/sub-XX/ses-YY/func      every SPM output goes here
%     D.anat    <preproc>/sub-XX/ses-YY/anat
%     D.fmap    <preproc>/sub-XX/ses-YY/fmap
%
%     D.roiSes       <rois>/sub-XX/ses-YY/anat
%     D.fieldmapSes  <fieldmaps>/sub-XX/ses-YY/fmap
%
%   WHY THIS EXISTS
%   ---------------
%   SPM writes next to its inputs, and worse, Realign & Unwarp and Coregister
%   rewrite the *headers of the input images themselves*. Run the pipeline
%   straight on a BIDS directory and the raw data quietly stops being raw.
%   So every file SPM is going to touch is copied here first, and the raw BIDS
%   tree is only ever read.
%
%   The three datasets are separate because they have different lifetimes:
%   spm-preproc is regenerated from scratch on every full preprocessing run
%   (see bids_reset_deriv), while rois and fieldmaps hold work that took a
%   human decision or a separate FSL run and must survive.
%
%   See also BIDS_RESET_DERIV, BIDS_STAGE.

if nargin < 2, SUB = ''; end
if nargin < 3, SES = ''; end

D = struct();

% ---- root ----
if isfield(cfgRaw, 'DERIV_ROOT') && ~isempty(cfgRaw.DERIV_ROOT)
    D.root = char(cfgRaw.DERIV_ROOT);
elseif isfield(cfgRaw, 'OUTPUT_ROOT') && ~isempty(cfgRaw.OUTPUT_ROOT)
    D.root = fullfile(char(cfgRaw.OUTPUT_ROOT), 'derivatives');
else
    error('bids_deriv:NoRoot', ...
        ['Cannot work out where the derivatives go.\n' ...
         'Set DERIV_ROOT, or OUTPUT_ROOT, in pipeline_config.cfg.']);
end

D.preproc   = fullfile(D.root, cfg_str(cfgRaw, 'DERIV_PREPROC',   'spm-preproc'));
D.rois      = fullfile(D.root, cfg_str(cfgRaw, 'DERIV_ROIS',      'rois'));
D.fieldmaps = fullfile(D.root, cfg_str(cfgRaw, 'DERIV_FIELDMAPS', 'fieldmaps'));

% ---- this session ----
D.sub = SUB;
D.ses = SES;

D.sesDir = session_path(D.preproc, SUB, SES);
D.func   = fullfile(D.sesDir, 'func');
D.anat   = fullfile(D.sesDir, 'anat');
D.fmap   = fullfile(D.sesDir, 'fmap');

D.roiSes      = fullfile(session_path(D.rois, SUB, SES), 'anat');
D.fieldmapSes = fullfile(session_path(D.fieldmaps, SUB, SES), 'fmap');
end

function p = session_path(base, SUB, SES)
p = base;
if ~isempty(SUB), p = fullfile(p, SUB); end
if ~isempty(SES), p = fullfile(p, SES); end
end

function v = cfg_str(cfgRaw, key, default)
if isfield(cfgRaw, key) && ~isempty(cfgRaw.(key))
    v = char(cfgRaw.(key));
else
    v = default;
end
end
