function data = eventDetection(x,y,ETparams)
%--------------------------------------------------------------------------
% README
%
% This code processes the recorded eye movement data to extract saccades,
% glissades and fixations.
% Observe that the algorithm is suitable ONLY for data collected from
% viewers keeping their heads relatively still while watching static
% stimuli. Note that at present it is not designed to handle data
% containing smooth pursuit movements.
%
% This code only requires only requires the statistics toolbox for the
% functions nanmean and nanstd. If you do not have access to this toolbox,
% you can get a nanmean implementation from Peter Acklam's website:
% http://home.online.no/~pjacklam/matlab/software/util/statutil/nanmean.m
% you'll have to implement nanstd yourself following this nanmean as an
% example

% Prepare data (move origin, things like that (TODO))
%-------------------------------------
data = prepareData(x,y,ETparams);

% Calculate velocity and acceleration
%-------------------------------------
data = filterDataAndCalcDerivative(data,ETparams);

% Detect blinks and noise
%-------------------------------------
data = removeNoise(data,ETparams);

% iteratively find the optimal noise threshold
%-------------------------------------
data = estimateSaccadeVelocityThresholds(data,ETparams);

% Detect saccades and glissades
% then, get information about them
%-------------------------------------
data = detectSaccadesAndGlissades (data,ETparams);
data = processSaccadesAndGlissades(data,ETparams);


% Implicitly detect fixations
%-------------------------------------
data = detectFixations (data,ETparams);
data = processFixations(data,ETparams);