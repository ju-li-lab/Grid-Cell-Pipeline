function C = spm_options()
% SPM_OPTIONS  Every MATLAB-side setting for the Grid Cell Pipeline.
%
% Paths, cluster resources and which subjects to process live in config.sh.
% Everything that ends up inside an SPM batch lives here.
%
% HOW THIS FILE WORKS
% -------------------
% The structs below are handed to SPM wholesale, e.g.
%
%     matlabbatch{1}.spm.spatial.realignunwarp.eoptions = C.realignunwarp.eoptions;
%
% so ANY option SPM accepts can be set by adding a field here. Nothing else
% needs changing - no parser, no mapping table, no validator to update.
% Field names are SPM's own; look them up in the SPM batch editor or in
% spm_cfg_realignunwarp.m / spm_cfg_fieldmap.m.
%
% Values are MATLAB literals, so [] is an empty number, '' is an empty
% string, and [12 12] is a vector. There is no guessing.
%
% One exception: the output filename prefixes (u, s, r) are NOT here. The
% shell scripts need them too, so config.sh owns them and the batch
% builders fill them in.

C = struct();


%% ========================================================================
%  1. ACQUISITION - must match your MRI protocol
%  ========================================================================
%  Wrong values here are silently wrong results. Check them against the
%  BIDS JSON sidecars in func/ and fmap/.

C.TR = 2.06;                    % Repetition time, SECONDS. From the task JSON.

% Echo times of the GRE fieldmap, MILLISECONDS.
% fmap/*_phasediff.json -> EchoTime1 and EchoTime2 (BIDS stores seconds).
TE_SHORT_MS = 4.92;
TE_LONG_MS  = 7.38;

% Total EPI readout, MILLISECONDS. Used by the SPM FieldMap toolbox (tert).
% func/*_bold.json -> TotalReadoutTime (seconds, so multiply by 1000).
%
% Readout time scales the distortion correction
% linearly, so one of the two is wrong by that factor. TotalReadoutTime and
% EstimatedTotalReadoutTime are different BIDS fields and may legitimately
% differ, but not usually by this much.
TOTAL_READOUT_MS = 33.2499;

% Blip direction for SPM FieldMap, +1 or -1. Depends on the phase-encode
% direction of the TASK, not the fieldmap.
% BIDS PhaseEncodingDirection: j -> typically +1, j- -> typically -1.
BLIP_DIRECTION = -1;

% Is the fieldmap itself EPI-based? 0 = no (GRE magnitude + phasediff),
% 1 = yes. A fieldmap made from a GRE phasediff is NOT EPI-based; with 1,
% SPM skips matching the fieldmap to the EPI.
EPI_BASED_FIELDMAP = 0;


%% ========================================================================
%  2. REALIGN & UNWARP  (spm.spatial.realignunwarp)
%  ========================================================================

C.realignunwarp.eoptions = struct( ...
    'quality', 1, ...
    'sep',     1.5, ...
    'fwhm',    1, ...
    'rtm',     1, ...
    'einterp', 4, ...
    'ewrap',   [0 1 0], ...
    'weight',  '');

C.realignunwarp.uweoptions = struct( ...
    'basfcn',   [12 12], ...
    'regorder', 1, ...
    'lambda',   100000, ...
    'jm',       0, ...          % Jacobian deformations: 0 = No, 1 = Yes.
    ...                         % Was 3, which is not a value SPM accepts.
    'fot',      [4 5], ...      % first-order effects: pitch and roll
    'sot',      [], ...         % second-order effects; [] = none
    'uwfwhm',   2, ...
    'rem',      1, ...
    'noi',      5, ...
    'expround', 'Average');

C.realignunwarp.uwroptions = struct( ...
    'uwwhich', [2 1], ...       % [2 1] = write resliced images + mean
    'rinterp', 6, ...
    'wrap',    [0 1 0], ...
    'mask',    0);              % prefix comes from config.sh (RESLICE_PREFIX)


%% ========================================================================
%  3. REALIGN ONLY  (spm.spatial.realign.estwrite)
%  ========================================================================
%  Used by PREPROC_MODE=realign_only, and after FSL applytopup.

C.realign.eoptions = struct( ...
    'quality', 1, ...
    'sep',     4, ...
    'fwhm',    5, ...
    'rtm',     1, ...
    'interp',  4, ...           % note: SPM calls this 'interp', not 'einterp'
    'wrap',    [0 1 0], ...
    'weight',  '');

C.realign.roptions = struct( ...
    'which',  [2 1], ...
    'interp', 6, ...
    'wrap',   [0 0 0], ...
    'mask',   0);               % prefix comes from config.sh (RESLICE_PREFIX)


%% ========================================================================
%  4. SPM FIELDMAP TOOLBOX  (spm.tools.fieldmap.calculatevdm)
%  ========================================================================
%  Builds the voxel displacement map in realign_unwarp and
%  precalc_fieldmap modes. Not used when PREPROC_MODE=topup.

C.fieldmap.defaultsval = struct( ...
    'et',        [TE_SHORT_MS TE_LONG_MS], ...
    'maskbrain', 1, ...
    'blipdir',   BLIP_DIRECTION, ...
    'tert',      TOTAL_READOUT_MS, ...
    'epifm',     EPI_BASED_FIELDMAP, ...
    'ajm',       0, ...         % apply Jacobian modulation
    ...
    ... % Phase unwrapping.
    'uflags', struct( ...
        'method', 'Mark3D', ...
        'fwhm',   5, ...
        'pad',    0, ...
        'ws',     1), ...
    ...
    ... % Brain mask SPM builds to unwrap the phase against. When GRE
    ... % unwarping looks wrong, this is the first thing to try: loosen
    ... % thresh, or change nerode/ndilate, and look at the resulting mask.
    'mflags', struct( ...
        'template', {{''}}, ...   % empty = SPM's own toolbox/FieldMap/T1.nii
        ...                      % (resolved by the VDM builder, once SPM is on the path)
        'fwhm',     5, ...
        'nerode',   2, ...
        'ndilate',  4, ...
        'thresh',   0.5, ...
        'reg',      0.02));

% How the VDM is written out.
C.fieldmap.vdmflags = struct( ...
    'udir',    1, ...
    'rinterp', 4, ...
    'wrap',    [0 1 0], ...
    'mask',    1, ...
    'order',   1, ...
    'prefix',  'vdm_');

% Match the VDM to the EPI it will be applied to (1) or not (0).
C.fieldmap.matchvdm      = 1;
C.fieldmap.writeunwarped = 0;
C.fieldmap.matchanat     = 0;


%% ========================================================================
%  5. FSL TOPUP  (PREPROC_MODE=topup)
%  ========================================================================

% Phase-encode direction of the BOLD task EPIs, FSL convention.
% BIDS PhaseEncodingDirection: i->x, i-->x-, j->y, j-->y-, k->z, k-->z-
C.topup.pe_dir_bold = 'y';

% ...and of the reverse-PE EPI, normally the opposite.
C.topup.pe_dir_reverse = 'y-';

% Total readout, SECONDS, for the acqparams file.
% func/*_bold.json -> EstimatedTotalReadoutTime. See the note on
% TOTAL_READOUT_MS above: these two currently disagree.
C.topup.readout_sec = 0.05409;

% Interpolation for applytopup: 'spline' or 'trilinear'
C.topup.interp = 'spline';


%% ========================================================================
%  6. COREGISTRATION  (spm.spatial.coreg.estwrite)
%  ========================================================================
%  T2w -> session mean EPI, then reslice the ROI masks into EPI space.

C.coreg.eoptions = struct( ...
    'cost_fun', 'nmi', ...
    'sep',      [4 2 1], ...
    'tol',      [0.02 0.02 0.02 0.001 0.001 0.001 ...
                 0.01 0.01 0.01 0.001 0.001 0.001], ...
    'fwhm',     [7 7]);

%  interp here applies to the ROI MASKS only:
%     0  nearest neighbour - one input sample per output voxel
%    -1  label (SPM25+)    - interpolates each label's indicator and takes
%                            the argmax, i.e. a local majority vote
%
%  KEEP THIS AT 0. Nearest neighbour is already very nearly volume-faithful
%  here: measured over all 36 masks in spm-preproc-oap, the true ROI volume
%  comes to 435 output voxels and NN produced 459 (+5%). A majority vote
%  comes to about 268 (-42%), because entorhinal cortex is a thin ribbon and
%  a voxel needs more than half its volume inside the ROI to survive the
%  vote. -1 was tried and gave smaller masks and fewer valid voxels.
%
%  The output grid is 2 x 2 x 3.45 mm (13.8 mm3), NOT 2 mm isotropic - 30
%  slices over 103.5 mm. That slice thickness, not the interpolation, is why
%  a ~150 mm3 ROI lands on only ~11 voxels.
%
%  The T2w is NOT resliced with this. coreg_reslice_rois.m writes it with
%  splines in a second call, because a T2w has a couple of thousand distinct
%  intensities and SPM's label path rejects it ("Too many label values").
C.coreg.roptions = struct( ...
    'interp', 0, ...
    'wrap',   [0 0 0], ...
    'mask',   0);               % prefix comes from config.sh (ROI_PREFIX)


%% ========================================================================
%  6b. SLICE TIMING  (spm.temporal.st)
%  ========================================================================
%  Where in the TR everything is shifted to, as a fraction: 0.5 is the
%  middle, 0 the first slice acquired.
%
%  This is tied to the GLM below. SPM samples the design matrix at
%  microtimeOnset of microtimeResolution bins, so after correcting to
%  fraction f you want microtimeOnset = f * microtimeResolution. At 0.5 and
%  16 bins that is 8, which is what section 9 has. The slicetime stage
%  checks the pair and warns if you change one without the other.
%
%  The slice TIMES themselves are not set here: they are read from each
%  run's JSON sidecar, so no slice order has to be assumed.
C.slicetime.ref_fraction = 0.5;


%% ========================================================================
%  7. SMOOTHING  (spm.spatial.smooth)
%  ========================================================================

% The kernel width and the prefix come from config.sh (SMOOTH_FWHM,
% SMOOTH_PREFIX), because the shell scripts need to know whether smoothing
% ran in order to find the right files afterwards.
C.smooth.dtype = 0;             % 0 = same datatype as input
C.smooth.im    = 0;             % no implicit masking


%% ========================================================================
%  8. SESSION MEAN
%  ========================================================================
%  The per-task means are averaged into one image, which is the
%  coregistration reference.

C.session_mean.enabled = true;
C.session_mean.name    = 'meanu_session.nii';


%% ========================================================================
%  9. GRIDCAT GLM
%  ========================================================================

C.GLM.xFoldSymmetry    = 6;     % 6 = hexagonal grid cells, 4 = square
C.GLM.maskingThreshold = 0.8;   % SPM masking threshold, 0-1. Higher = stricter.
C.GLM.microtimeOnset     = 8;
C.GLM.microtimeResolution = 16;
C.GLM.derivatives      = [0 0]; % HRF derivatives: [0 0] none,
                                % [1 0] temporal, [1 1] temporal + dispersion
C.GLM.HPF_cutoff       = 128;   % high-pass filter, seconds
C.GLM.dispDesignMatrix = 0;     % 1 to show the design matrix (debugging only)

C.GLM.eventUsageSpecifier_GLM1 = 2;
C.GLM.eventUsageSpecifier_GLM2 = 3;
C.GLM.keepUnusedGridEvents     = 1;
C.GLM.useVoxelWeighting        = 1;
C.GLM.avgOriAcrossRuns         = 0;

% Which regressor to use for grid events:
%   'pmod'                        one regressor with a parametric modulation
%   'aligned_misaligned'          one regressor for events aligned with the
%                                 mean grid orientation, one for misaligned
%   'aligned_misaligned_multiple' one regressor per orientation
C.GLM.gridRegressorMethod = 'pmod';

end
