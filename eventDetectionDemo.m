clear all, close all; clc
%%-------------------------------------------------------------------------
%%% This code is an implementation of 
%%% Nyström, M. & Holmqvist, K. (in press), "An adaptive algorithm for
%%% fixation, saccade, and glissade detection in eye-tracking data".
%%% Behavior Research Methods
%%-------------------------------------------------------------------------
addpath(genpath(fullfile(cd,'utils')))
addpath(genpath(fullfile(cd,'utilsExternal')))
%%-------------------------------------------------------------------------
%%% Init parameters
%%-------------------------------------------------------------------------

% load parameters
ETparams = getParameters;

% process params
ETparams = prepareParameters(ETparams);


%%-------------------------------------------------------------------------
%%% Begin detection
%%-------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Process eye-movement data, per trial
%--------------------------------------------------------------------------
[files,nfiles] = FileFromFolder(fullfile(cd,'example data'),'','txt');
fhndl = -1;
for p = 1:nfiles
    % Process data
    disp(files(p).fname)
    data = readNumericFile(fullfile(cd,'example data',files(p).name),4);
    results = eventDetection(data(:,2),data(:,3),data(:,4),ETparams,data(:,1));
    
    if 1
        % plot the trial (eye X, eye Y, velocity traces and scanpath,
        % as well as detected events)
        if ~ishghandle(fhndl)
            fhndl = figure('Units','normalized','Position',[0 0 1 1]);  % make fullscreen figure
        else
            clf(fhndl);
        end
        plotDetection(results,'deg','vel',ETparams.samplingFreq,ETparams.glissade.searchWindow,ETparams.screen.rect,'title',files(p).fname);
        set(fhndl,'Visible','on');  % assert visibility to bring window to front again after keypress
        pause
        if ~ishghandle(fhndl)
            break;
        end
    end
end
