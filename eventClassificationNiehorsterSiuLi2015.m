clear all, close all; clc
%%-------------------------------------------------------------------------
%%% This code is an implementation of 
%%% Nyström, M. & Holmqvist, K. (2010), "An adaptive algorithm for
%%% fixation, saccade, and glissade detection in eye-tracking data".
%%% Behavior Research Methods 42(1):188-204.
%%% https://doi.org/10.3758/BRM.42.1.188
%%%
%%% Includes multiple extensions developed by DN that were not in the
%%% article.
%%-------------------------------------------------------------------------
addpath(genpath(fullfile(cd,'function_library')))
addpath(genpath(fullfile(cd,'post-process')))       % post-process folder needed for Niehorster, Siu & Li (2015) code
%%-------------------------------------------------------------------------
%%% Init parameters
%%-------------------------------------------------------------------------

% settings for this file
plotWhat        = 'deg';                            % 'deg' or 'pix': determines in what unit to plot data
whichVel        = 'velX';                           % 'vel', 'velX' or 'velY': determines which velocity trace to plot when interface opens (there are buttons to view the others)
doPlot          = true;

% load parameters for event classifier
ETparams = defaultParameters;

% change some defaults as needed for this analysis:
ETparams.data.alsoStoreComponentDerivs  = true;
ETparams.data.detrendWithMedianFilter   = true;
ETparams.data.applySaccadeTemplate      = true;
ETparams.data.minDur                    = 500;
ETparams.fixation.doClassify            = false;
ETparams.blink.replaceWithInterp        = true;
ETparams.blink.replaceVelWithNan        = true;

% some more just to demo functionality
ETparams.saccade.useTemplateForRefine   = false;    % default false: if true, trace produced by cross-correlation of 2D velocity trace with saccade template is used for saccade onset and offset refinement
if strcmp(plotWhat,'pix')
    ETparams.data.alsoStoreandDiffPixels    = true;     % if true, velocity and acceleration in pixels is also stored, so can be plotted
end

% settings for code specific to Niehorster, Siu & Li (2015)
extraCut    = [0 50];                       % extra ms of data to cut before and after saccade.
qInterpMissingPos   = true;                 % interpolate using straight lines to replace missing position signals?

% settings for the saccade cutting (see cutSaccades.m for documentation)
cutPosTraceMode     = 2;
cutVelTraceMode     = 1;
cutSaccadeSkipWindow= 1;    % don't cut during first x seconds

% process params
ETparams = prepareParameters(ETparams);


%%-------------------------------------------------------------------------
%%% Begin classification
%%-------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Process eye-movement data, per trial
%--------------------------------------------------------------------------
[files,nfiles] = FileFromFolder(fullfile(cd,'example data','NiehorsterSiuLi2015'),'','txt');
fhndl = -1;
for p = 1:nfiles
    %%% 1. event classification and standard plot
    % read in data
    disp(files(p).fname)
    rawData = readNumericFile(fullfile(cd,'example data','NiehorsterSiuLi2015',files(p).name),4);
    
    % perform event classification
    data    = runNH2010Classification(rawData(:,2),rawData(:,3),rawData(:,4),ETparams,rawData(:,1));
    
    %%% 2. for Niehorster, Siu & Li (2015), we had to remove saccades from
    %%% the data traces, do that here and plot again.
    % merge all glissades with the saccades (timestamp manipulation...)
    data = mergeSaccadesAndGlissades(data);
    if isfield(data,'glissade')
        data = rmfield(data,'glissade');    % remove glissade information, we merged them anyway
    end
    
    % replace missing data by linearly interpolating position and velocity
    % between start and end of each missing interval (so, creating a ramp
    % between start and end position/velocity).
    data = replaceMissing(data,qInterpMissingPos);
    
    % desaccade velocity and/or position
    data = cutSaccades(data,ETparams,cutPosTraceMode,cutVelTraceMode,extraCut,cutSaccadeSkipWindow);
    % construct saccade only traces
    data = cutPursuit(data,ETparams,1);
    
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
