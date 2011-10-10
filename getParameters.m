function ETparams = getParameters

% user settings
ETparams.screen.resolution              = [1280 1024];
ETparams.screen.size                    = [0.40 0.30];
ETparams.screen.viewingDist             = 0.58;
ETparams.screen.subjectStraightAhead    = [640 381.5];  % Specify the screen coordinate that is straight ahead of the subject. Just specify the middle of the screen unless its important to you to get this very accurate!

% flip the Y coordinate of the data? All the routines assume the origin of
% the screen (0,0) is at the top left corner. You'll have to flip if the
% your data's origin is the lower left corner. Do a flip X if your origin
% is on the right side of the screen (sic).
ETparams.data.qFlipY                    = false;
ETparams.data.qFlipX                    = false;
% By default, eye velocity and acceleration are computed with a
% Savitzky-Golay differentiating filter (it basically fits a second order
% polynomial in a moving window to the eye position data and takes the
% derivatives analytically). Set the below to true to simply numerical
% differentiate with matlab's diff()
ETparams.data.qNumericallyDifferentiate = true;
% If true, datastream of eyeposition in pixels is also stored, smoothed and
% the derivatives are taken. Might be needed in some usage cases. The
% eventDetection however always runs on eye position in degrees.
ETparams.data.qAlsoStoreandSmoothPixels = true;
% Do a precise calculation of angular eye velocity and acceleration? If
% not, we compute derivatives of eye azimuth and elevation analytically
% from the parameters of a fitted polynomial and then apply Pythagoras'
% theorem to compute eye velocity/acceleration. This is crude and should
% not be used if you're interested in the eye velocity, but its sufficient
% if you simply want to detect saccades in periods of fixation and/or
% smooth pursuit and are not interested in accurate measures of eye
% velocity/acceleration.
ETparams.data.qPreciseCalcDeriv         = false;
ETparams.data.qAlsoStoreComponentDerivs = true;         % if true, velocity in X/azimuth and Y/elevation direction separately are also stored.

% Option to use median filter for detrending velocity data (e.g. removing
% pursuit baseline speed). This is only useful if saccade templates are
% used, detrended velocity is then used as input to xcorr with the saccade
% template
ETparams.data.qDetrendWithMedianFilter  = true;
ETparams.data.qDetrendAll               = false;        % if true, all velocity traces (also pixels and also components, if available) will be detrended. If false, only eyevelocity in degrees will be done. Only set this to true if you want this data, this code doesn't use it
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
ETparams.data.qApplySaccadeTemplate     = true;

ETparams.samplingFreq                   = 500;

ETparams.blink.velocityThreshold        = 1000;         % if vel > 1000 °/s, it is noise or blinks
ETparams.blink.accThreshold             = 100000;       % if acc > 100000 °/s², it is noise or blinks

ETparams.saccade.peakVelocityThreshold  = 100;          % Initial value of the peak detection threshold, °/s
ETparams.saccade.xCorrPeakThreshold     = .2;           % Initial threshold for saccade detection from data filtered by saccade template
ETparams.saccade.qSaccadeTemplateRefine = false;        % saccade beginnings and ends are refined from the xcorr response of the saccade template, not from the velocity trace
ETparams.saccade.minDur                 = 10;           % in milliseconds
ETparams.saccade.allowNaN               = true;         % if true, allow NaNs in saccade intervals

ETparams.glissade.qDetect               = true;         % if true, do glissade detection
ETparams.glissade.searchWindow          = 40;           % window after saccade in which we search for glissades, in milliseconds
ETparams.glissade.maxDur                = 80;           % in milliseconds
ETparams.glissade.allowNaN              = false;        % if true, allow NaNs in saccade intervals

% fixation here is defined as 'not saccade or glissade', it could thus also
% be pursuit
ETparams.fixation.qDetect               = false;        % if true, do fixation detection
ETparams.fixation.minDur                = 100;          % in milliseconds (long as we're not interested in very short pursuit intervals)
% How to deal with NaNs during possible fixation periods:
% 1: do not allow NaN during fixations, whole fixation thrown out
% 2: ignore NaNs and calculate mean fixation position based on other data
%    (not recommended in almost any situation, if you don't like 1,
%    consider option 3)
% 3: split fixation into multiple, providing each is at least minDur long
%    (e.g. one 250 ms fixation with some data missing in the middle might
%    be split up into a 100 ms and a 120 ms fixation)
ETparams.fixation.treatNaN              = 1;