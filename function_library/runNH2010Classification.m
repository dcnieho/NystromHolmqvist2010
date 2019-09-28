function [data,ETparams] = runNH2010Classification(x,y,pupilsize,ETparams,varargin)
%--------------------------------------------------------------------------
% README
%
% This code is based on an implementation of Nyström, M. & Holmqvist, K.
% (2010), "An adaptive algorithm for fixation, saccade, and glissade
% detection in eye-tracking data". Behavior Research Methods 42(1):188-204.
% It processes the recorded eye movement data to extract saccades,
% fixations, and glissades (the latter are now recognized to be
% post-saccadic oscillations).
%
% Please note that there have been several extensions to this algorithm to
% have it optionally allow some missing data during events, work better
% during pursuit, do blink classification, etc.

% Prepare data and params (move origin, things like that)
%-------------------------------------
[data,ETparams] = prepareDataAndParams(x,y,pupilsize,ETparams,varargin{:});

% Calculate velocity and acceleration
%-------------------------------------
data = filterDataAndCalcDerivative(data,ETparams);

% Classify noise episodes
%-------------------------------------
data = removeNoise(data,ETparams);

% Detrend and apply saccade template
%-------------------------------------
data = detrendAndApplySaccadeTemplate(data,ETparams);

% Iteratively find the optimal classification threshold
%-------------------------------------
data = estimateThresholds(data,ETparams);

% Classify saccades and glissades
%-------------------------------------
data = classifySaccadesAndGlissades(data,ETparams);

% Classify and remove blinks
% Also classifies some noisy bits of data as evidenced
% by wobbly pupil size trace. Thats ok, err on save
% side here, make sure we don't miss any even small
% eyelid droops, and some noise removal is nice too...
% (Need to do this before saccades are merged, or
% we might also remove saccades that occured right
% after a blink)
%-------------------------------------
if ETparams.blink.classifyMode
    data = classifyAndRemoveBlinks(data,ETparams);
end

% Identify episodes of missing data in the traces
%-------------------------------------
data = flagMissing(data);

% Now merge saccades with short intervals between them
% and get information about them
%-------------------------------------
data = processSaccadesAndGlissades(data,ETparams);

% Implicitly classify fixations
% and get information about them
%-------------------------------------
if ETparams.fixation.doClassify
    data = classifyFixations (data,ETparams);
    data = processFixations(data,ETparams);
end
