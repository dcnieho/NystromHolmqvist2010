function ETparams = getParameters

% user settings
ETparams.screen.resolution              = [1680 1050];
ETparams.screen.size                    = [0.473 0.296];
ETparams.screen.viewingDist             = 0.60;
ETparams.screen.dataCenter              = [840 525];    % center of screen has these coordinates in data
ETparams.screen.subjectStraightAhead    = [840 525];    % Specify the screen coordinate that is straight ahead of the subject. Just specify the middle of the screen unless its important to you to get this very accurate!

% flip the Y coordinate of the data? All the routines assume the origin of
% the screen (0,0) is at the top left corner. You'll have to flip if the
% your data's origin is the lower left corner. Do a flip X if your origin
% is on the right side of the screen (sic).
ETparams.data.qFlipY                    = true;
ETparams.data.qFlipX                    = false;
% By default, eye velocity and acceleration are computed with a
% Savitzky-Golay differentiating filter (it basically fits a second order
% polynomial in a moving window to the eye position data and takes the
% derivatives analytically). Set the below qNumericallyDifferentiate to
% true to do a simple numerical difference by convoluting with [1 0 -1] (a
% slightly smoothed and correctly centered version of matlab's diff).
ETparams.data.qNumericallyDifferentiate = true;         % if 1, do above convolution methods (moveing average of three samples), if 2, do simple diff()
ETparams.data.filterWindow              = 10;           % ms, if using Savitzky-Golay, filter window length. Make sure its narrower than smallest features you want to detect
ETparams.data.qAlsoStoreComponentDerivs = true;         % if true, velocity in X/azimuth and Y/elevation direction separately are also stored.
% If true, eyeposition trace in pixels is also stored and derivatives
% (smoothed, if using Savitzky-Golay) are calculated. Might be needed in
% some usage cases. The eventDetection however always runs on eye
% pos/vel/acc in degrees.
ETparams.data.qAlsoStoreandDiffPixels   = true;

% Option to use median filter for detrending velocity data (e.g. removing
% pursuit baseline speed). This is only useful if saccade templates are
% used. Detrended velocity, if available, is then used as input to xcorr
% with the saccade template.
ETparams.data.qDetrendWithMedianFilter  = false;
ETparams.data.qDetrendAll               = false;        % if true, all velocity traces (also pixels and also components, if available) will be detrended. If false, only 2D eye velocity in degrees will be done. Only set this to true if you want this data, the code doesn't use it
ETparams.data.medianWindowLength        = 40;           % ms

% Option to first convolve velocity trace with the velocity profile of a
% saccade. This will bring out the saccades by reducing the amplitude of
% features in the trace that are not like saccades. If this is set to true,
% saccadic peaks are identified based on this xcorr response. Whether
% refinement of saccade starts and ends, and glissade detection is also
% done on the response trace however depends on
% ETparams.saccade.qSaccadeTemplateRefine. I'd recommend to leave that to
% false as the profiles of the saccades is distorted after convolution with
% the template.
% See also ETparams.saccade.xCorrPeakThreshold
ETparams.data.qApplySaccadeTemplate     = false;

ETparams.data.minDur                    = 1000/60;      % in milliseconds, minimum stretch of consequtive data points. If a shorter stretch is found in between missing data, it is flagged as missing as well.
ETparams.data.missCoords                = [0 0];        % if empty, not used, else if tracking indicates missing by sending specific coordinates, put them here

ETparams.samplingFreq                   = 120;

% blink detection. Note that there are other ways to detect blinks in the
% code, e.g., all saccades are checked for their level of blinkness
ETparams.blink.detectMode               = 1;            % if >0, do blink detection. 1: only use thresholding of pupil size change. 2: only use vel/acc thresholding. 3: do both to detect potential blinks
ETparams.blink.dSizeThreshold           = 25000;        % Initial threshold for blink detection from pupil size change data
ETparams.blink.nStd                     = 5;            % threshold is set at mean pupil size change + nStd*std of pupil size change by optimization algorithm
ETparams.blink.localNoiseWindowLength   = 50;           % in milliseconds, window before a blink in which to calculate noise and mean pupil size change, used to calculate blink offset thresholds
ETparams.blink.minPeakSamples           = 2;            % minimum number of samples data need to be above peak threshold for a peak to be considered a potential blink. Very short peaks are likely to be noise.
ETparams.blink.qExcludeOneSampleBlinks  = true;         % Based on pupil size (instead of change of size): don't consider single sample eye closed as a blink (might be multiple single asmaple occurences)
ETparams.blink.minDur                   = 80;           % in milliseconds
ETparams.blink.velocityThreshold        = 1000;         % if vel > 1000 °/s, it is noise or blinks
ETparams.blink.accThreshold             = 100000;       % if acc > 100000 °/s², it is noise or blinks
ETparams.blink.mergeWindow              = 75;           % merge blinks that are less than this many non-NaN samples apart. Set to 0 if you don't want any merging.
ETparams.blink.qGrowToOverlapSaccade    = true;         % if blink detected, see if detected on and offsets overlap with saccades. If so, take on-/offsets of these saccades as the blink on-/offsets
ETparams.blink.qGrowToOverlapGlissade   = false;        % if blink detected, if blink.qGrowToOverlapSaccade is true, see if overlapping saccade has a glissade, and if so grow to end of glissade. If blink.qGrowToOverlapSaccade is false, see if detected on and offsets overlap with glissades. If so, take on-/offsets of these glissades as the blink on-/offsets
ETparams.blink.qReplaceWithInterp       = false;        % replace position, velocity, acceleration with linear interpolation between bounds. overrides ETparams.blink.qReplaceVelWithNan
ETparams.blink.qReplaceVelWithNan       = false;        % if true, blinks in all the velocity traces are replaced with NaN

ETparams.saccade.peakVelocityThreshold  = 100;          % Initial value of the peak detection threshold, °/s
ETparams.saccade.peakXCorrThreshold     = .2;           % Initial threshold for saccade detection from data filtered by saccade template
ETparams.saccade.qSaccadeTemplateRefine = false;        % saccade beginnings and ends are refined from the xcorr response of the saccade template (true), not from the velocity trace (false). Leave this to false to avoid distortion of saccade beginning and end due to low-pass filtering of template. Must be false if saccade template isn't used
ETparams.saccade.localNoiseWindowLength = 50;           % in milliseconds, window before a saccade in which to calculate noise and mean eye speed, used to calculate saccade offset thresholds
ETparams.saccade.minPeakSamples         = 1;            % minimum number of samples data need to be above peak threshold for a peak to be considered a potential saccade. Very short peaks are likely to be noise.
ETparams.saccade.minDur                 = 10;           % in milliseconds
ETparams.saccade.mergeWindow            = -1;           % merge saccades that are less than this apart (this is counted from saccade or glissade end (if any) to next saccade start). Set to 0 if you don't want any merging.
ETparams.saccade.seWindowSamp           = 4;            % number of samples before onset and after offset to use for calculating saccade start and end points (onset and offset themselves are always used, this is the number of extra samples)
ETparams.saccade.allowNaN               = true;         % if true, allow NaNs in saccade intervals. If false, blink detection dies as its basic assumption is that blinks were already detected as saccades due to their large vertical velocity

ETparams.glissade.qDetect               = true;         % if true, do glissade detection
ETparams.glissade.searchWindow          = 40;           % window after saccade in which we search for glissades, in milliseconds
ETparams.glissade.maxDur                = 80;           % in milliseconds
ETparams.glissade.seWindowSamp          = 0;            % number of samples before onset and after offset to use for calculating glissade start and end points (onset and offset themselves are always used, this is the number of extra samples)
ETparams.glissade.allowNaN              = false;        % if true, allow NaNs in saccade intervals

% fixation here is defined as 'not saccade or glissade', it could thus also
% be pursuit. For detecting "fixations", the code starts out with all
% intervals that are not saccade or glissade and then processes them
% according to the below settings. It could thus be that some parts of the
% trace end up unclassified as they don't qualify as saccade, glissade or
% "fixation".
ETparams.fixation.qDetect               = false;        % if true, do fixation detection
ETparams.fixation.minDur                = 70;           % in milliseconds
% How to deal with NaNs during possible fixation periods:
% 1: do not allow NaN during fixations, whole fixation thrown out
% 2: ignore NaNs and calculate mean fixation position based on available
%    data. Checks that no significant eye position jumps occur during
%    missing data
% 3: split fixation into multiple as delimited by nan, providing each is at
%    least fixation.minDur long (e.g. one 250 ms fixation with some data
%    missing in the middle might be split up into a 100 ms and a 120 ms
%    fixation)
ETparams.fixation.treatNaN              = 1;
ETparams.fixation.NaNMaxJump            = 1;            % maximum amplitude (degree) of eye fixation position change allowed during missing data (if ETparams.fixation.treatNaN==2)



%%%%%%%%%%%% Ignore / leave to false for now:
% Do a precise calculation of angular eye velocity and acceleration? If
% not, we apply Pythagoras' theorem to compute eye velocity/acceleration
% from the azimuthal and elevational coordinate velocities. When we have no
% knowledge of torsional eye movements and therefore have to assume that
% they are 0, both methods are equivalent. The precise calculations then
% still give you the axis of rotation as well, but as we assume torsion is
% 0, this axis is systematically biased for all but the smallest movements
% away from the primary reference position. So there is no use in using the
% precise calculations when you don't have eye torsion information. So when
% you are interested in exact eye velocities, you'd do well to acquire
% measures of torsion as well. Nonetheless, this straightforward
% calculation of 2D eye velocity is sufficient for accurate detection of
% saccades if that is all that you are interested in.
ETparams.data.qPreciseCalcDeriv         = false;
