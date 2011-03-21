function ETparams = getParameters

% user settings
ETparams.screen.resolution              = [1280 1024];
ETparams.screen.size                    = [0.38 0.30];
ETparams.screen.viewingDist             = 0.67;
ETparams.screen.subjectStraightAhead    = [640 312];    % Specify the screen coordinate that is straight ahead of the subject. Just specify the middle of the screen unless its important to you to get this very accurate!

% flip the Y coordinate of the data? All the routines assume the origin of
% the screen (0,0) is at the top left corner. You'll have to flip if the
% your data's origin is the lower left corner. Do a flip X if your origin
% is on the right side of the screen (sic).
ETparams.data.qFlipY                    = false;
ETparams.data.qFlipX                    = false;
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
ETparams.data.qAlsoStoreComponentDerivs = true;         % if true, velocity in X and Y direction separately are also stored.

ETparams.samplingFreq                   = 500;

ETparams.blink.velocityThreshold        = 1000;         % if vel > 1000 °/s, it is noise or blinks
ETparams.blink.accThreshold             = 100000;       % if acc > 100000 °/s², it is noise or blinks

ETparams.saccade.peakVelocityThreshold  = 100;          % Initial value of the peak detection threshold, °/s
ETparams.saccade.minDur                 = 10;           % in milliseconds
ETparams.saccade.allowNaN               = true;         % if true, allow NaNs in saccade intervals

ETparams.glissade.qDetect               = true;         % if true, do glissade detection
ETparams.glissade.searchWindow          = 40;           % window after saccade in which we search for glissades, in milliseconds
ETparams.glissade.maxDur                = 80;           % in milliseconds
ETparams.glissade.allowNaN              = false;        % if true, allow NaNs in saccade intervals

ETparams.fixation.qDetect               = false;        % if true, do fixation detection
ETparams.fixation.minDur                = 40;           % in milliseconds
% How to deal with NaNs during possible fixation periods:
% 1: do not allow NaN during fixations
% 2: ignore NaNs and calculate mean fixation position based on other data
%    (not recommended in almost any situation, if you don't like 1,
%    consider option 3)
% 3: split fixation into multiple, providing each is at least minDur long
%    (e.g. one 250 ms fixation with some data missing in the middle might
%    be split up into a 100 ms and a 120 ms fixation)
ETparams.fixation.treatNaN              = 1;