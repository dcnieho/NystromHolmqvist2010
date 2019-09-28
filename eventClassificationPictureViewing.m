clear all, close all; clc
%%-------------------------------------------------------------------------
%%% This code is an implementation of Nyström, M. & Holmqvist, K. (2010).
%%% It includes multiple extensions developed by DN that are listed in the
%%% readme.md file. When using this code, please cite Niehorster, Siu & Li
%%% (2015). If using ETparams.saccade.onsetRefineMethod=2, please
%%% additionally cite Oliva, Niehorster, Jarodzka & Holmqvist (2017).
%%%
%%% Example citation:
%%% Saccades were classified using the Niehorster, Siu & Li (2015)
%%% implementation of the Nyström & Holmqvist (2010) algorithm, with
%%% default settings. In addition, saccade onsets were determined using the
%%% method of Oliva, Niehorster, Jarodzka & Holmqvist (2017).
%%%
%%% References:
%%% Nyström, M. & Holmqvist, K. (2010), "An adaptive algorithm for
%%%    fixation, saccade, and glissade detection in eye-tracking data".
%%%    Behavior Research Methods 42(1): 188-204. doi: 10.3758/BRM.42.1.188
%%% Niehorster, D.C., Siu, W.W.F., & Li, L. (2015). Manual tracking
%%%    enhances smooth pursuit eye movements. Journal of Vision 15(15), 11.
%%%    doi: 10.1167/15.15.11
%%% Oliva, M., Niehorster, D.C., Jarodzka, H., & Holmqvist, K. (2017).
%%%    Social Presence Influences Saccadic and Manual Responses.
%%%    I-Perception 8(1). doi: 10.1177/2041669517692814
%%-------------------------------------------------------------------------
addpath(genpath(fullfile(cd,'function_library')))
%%-------------------------------------------------------------------------
%%% Init parameters
%%-------------------------------------------------------------------------

% settings for this file
plotWhat        = 'deg';                            % 'deg' or 'pix': determines in what unit to plot data
whichVel        = 'vel';                            % 'vel', 'velX' or 'velY': determines which velocity trace to plot when interface opens (there are buttons to view the others)
doPlot          = true;

% load parameters for event classifier
ETparams = defaultParameters;

% change some defaults as needed for this analysis:
ETparams.samplingFreq                   = 1000;
ETparams.screen.resolution              = [1920 1080];
ETparams.screen.size                    = [0.53 0.30];
ETparams.screen.viewingDist             = 0.80;
ETparams.screen.dataCenter              = [ 960  540];
ETparams.screen.subjectStraightAhead    = [ 960  540];

% some more just to demo functionality
ETparams.data.alsoStoreComponentDerivs  = true;         % set to true so we have 1D velocities and accelerations in the output as well
ETparams.data.applySaccadeTemplate      = false;        % default false, set to true to run event classification on a trace produced by the cross-correlation of gaze velocity with a saccade template, instead of on gaze velocity directly
if strcmp(plotWhat,'pix')
    ETparams.data.alsoStoreandDiffPixels    = true;     % if true, velocity and acceleration in pixels is also stored, so can be plotted
end

% process params
ETparams = prepareParameters(ETparams);


%%-------------------------------------------------------------------------
%%% Begin classification
%%-------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Process eye-movement data, per trial
%--------------------------------------------------------------------------
[files,nfiles] = FileFromFolder(fullfile(cd,'example data','pictureViewing'),'','txt');
fhndl = -1;
for p = 1:nfiles
    %%% 1. event classification and standard plot
    % read in data
    disp(files(p).fname)
    rawData = readNumericFile(fullfile(cd,'example data','pictureViewing',files(p).name),7,1);
    
    % perform event classification (arbitrarily on left eye
    data    = runNH2010Classification(rawData(:,2),rawData(:,3),rawData(:,6),ETparams,rawData(:,1));
    
    if doPlot
        % plot the trial (eye X, eye Y, velocity traces and scanpath,
        % as well as classified saccades and fixations)
        if ~ishghandle(fhndl)
            fhndl = figure('Units','normalized','Position',[0 0 1 1]);  % make fullscreen figure
        else
            clf(fhndl);
        end
        plotClassification(data,plotWhat,whichVel,ETparams.samplingFreq,ETparams.glissade.searchWindow,ETparams.screen.rect,'title',files(p).fname,'zeroStartT',true);
        set(fhndl,'Visible','on');  % assert visibility to bring window to front again after keypress
        pause
        if ~ishghandle(fhndl)
            break;
        end
    end
end

% clean up
if ishghandle(fhndl)
    close(fhndl);
end
rmpath(genpath(fullfile(cd,'function_library')))