function [data,ETparams] = runNH2010Classification(x,y,pupilsize,ETparams,varargin)
%--------------------------------------------------------------------------
% README
%
% This code is an implementation of Nyström, M. & Holmqvist, K. (2010).
% It includes multiple extensions developed by DN that are listed in the
% readme.md file. When using this code, please cite Niehorster, Siu & Li
% (2015). If using ETparams.saccade.onsetRefineMethod=2, please
% additionally cite Oliva, Niehorster, Jarodzka & Holmqvist (2017).
%
% Example citation:
% Saccades were classified using the Niehorster, Siu & Li (2015)
% implementation of the Nyström & Holmqvist (2010) algorithm, with
% default settings. In addition, saccade onsets were determined using the
% method of Oliva, Niehorster, Jarodzka & Holmqvist (2017).
%
% References:
% Nyström, M. & Holmqvist, K. (2010), "An adaptive algorithm for
%    fixation, saccade, and glissade detection in eye-tracking data".
%    Behavior Research Methods 42(1): 188-204. doi: 10.3758/BRM.42.1.188
% Niehorster, D.C., Siu, W.W.F., & Li, L. (2015). Manual tracking
%    enhances smooth pursuit eye movements. Journal of Vision 15(15), 11.
%    doi: 10.1167/15.15.11
% Oliva, M., Niehorster, D.C., Jarodzka, H., & Holmqvist, K. (2017).
%    Social Presence Influences Saccadic and Manual Responses.
%    I-Perception 8(1). doi: 10.1177/2041669517692814
%--------------------------------------------------------------------------

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
