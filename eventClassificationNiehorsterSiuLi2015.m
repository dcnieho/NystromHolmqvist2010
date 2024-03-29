clear all, close all; clc
%%-------------------------------------------------------------------------
%%% This code is an implementation of Nystr�m, M. & Holmqvist, K. (2010).
%%% It includes multiple extensions developed by DN that are listed in the
%%% readme.md file. When using this code, please cite Niehorster, Siu & Li
%%% (2015). If using ETparams.saccade.onsetRefineMethod=2, please
%%% additionally cite Oliva, Niehorster, Jarodzka & Holmqvist (2017).
%%%
%%% Example citation:
%%% Saccades were classified using the Niehorster, Siu & Li (2015)
%%% implementation of the Nystr�m & Holmqvist (2010) algorithm, with
%%% default settings. In addition, saccade onsets were determined using the
%%% method of Oliva, Niehorster, Jarodzka & Holmqvist (2017).
%%%
%%% References:
%%% Nystr�m, M. & Holmqvist, K. (2010), "An adaptive algorithm for
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
            % make fullscreen figure
            fhndl = figure('Units','pixels');
            if isprop(fhndl,'WindowState')
                fhndl.WindowState = 'Maximized';
                drawnow
            else
                ws    = get(0,'ScreenSize');
                hmmar = [0 0 0 40];    % left right top bottom -- works well for Windows 10
                fhndl.OuterPosition = [ws(1) + hmmar(1), ws(2) + hmmar(4), ws(3)-hmmar(1)-hmmar(2), ws(4)-hmmar(3)-hmmar(4)];
            end
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
rmpath(genpath(fullfile(cd,'post-process')))