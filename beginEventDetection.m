clear all, clear functions, close all; clc

% TODO notes
% - center data so that middle of screen is (0,0)
% - add extra pass to noise deletion, deleting sections of data that are
%   too short to be useful (< minfixdur)
% - name all function using detectXX, processXX and ComputeXX
% - split up in multiple parts: eyemovement detection is separate from
%   analysis, saved to file in between (as old code)
% -   therefore do not have to worry about the current setup with 
% - add switch for calculating precise eye position, based on cosine rule
% - port functionality of filternanfix to here, but allow to switch it on
%   or off by a boolean
% - add setting for how many samples sign of derivative of velocity needs
%   to change to signal offset (i suggest 3, but test carefully)

%%-------------------------------------------------------------------------
%%% Init parameters
%%-------------------------------------------------------------------------
load('1250Hz_3_Participants.mat');

% user settings
ETparams.screenSz = [1024 768];
ETparams.screenDim = [0.38 0.30];
ETparams.viewingDist = 0.67;
ETparams.samplingFreq = 1250;
ETparams.blinkVelocityThreshold = 1000;             % if vel > 1000 degrees/s, it is noise or blinks
ETparams.blinkAccThreshold = 100000;                % if acc > 100000 degrees/s^2, it is noise or blinks
ETparams.peakDetectionThreshold = 100;              % Initial value of the peak detection threshold.

ETparams.minSaccadeDur = 10;                        % in milliseconds
ETparams.glissadeSearchWindow = 40;                 % window after saccade in which we search for glissades, in milliseconds
ETparams.maxGlissadeDur = 80;                       % in milliseconds
ETparams.minFixDur = 40;                            % in milliseconds

% process them
% Calculate how many degrees one pixel spans.
[ETparams.angleInPixelsH, ETparams.angleInPixelsV] = degrees2pixels(1, ETparams.viewingDist, ETparams.screenSz, ETparams.screenDim);

%%-------------------------------------------------------------------------
%%% Begin detection
%%-------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Process eye-movement data for file (participant) and trial separately, i - files, j -
% trials
%--------------------------------------------------------------------------
for i = 2%1:size(ETparams.data,1)
    for j = 1%1:size(ETparams.data,2)
        % Process data
        fprintf('Subj %d, Trial %d\n',i,j);
        data = eventDetection(ETdata(i,j).X,ETdata(i,j).Y,ETparams);
        
        if 1
            % plot the trial (eye X, eye Y, velocity traces and scanpath,
            % as well as detected events (TODO: scanpath)
            data.sampleRate = ETparams.samplingFreq;
            data.glissadeSearchWindow = ETparams.glissadeSearchWindow;
            fhndl = figure('Units','normalized','Position',[0 0 1 1]);  % make fullscreen figure
            plotDetection(data);
            set(fhndl,'Visible','on');  % assert visibility to bring window to front again after keypress
            pause
            if ~ishghandle(fhndl)
                break;
            end
        end
    end
end

save([cd,'\DetectionResults\DetectionResults.mat'],'ETparams');

%%-------------------------------------------------------------------------
%%% Plot results
%%-------------------------------------------------------------------------

% Calculate basic parameters
mean(cat(1,ETparams.glissadeInfo.duration))
mean(cat(1,ETparams.saccadeInfo.duration))
mean(cat(1,ETparams.fixationInfo.duration))

%Plot histograms
% figure
% hist(cat(1,ETparams.glissadeInfo.duration),100)
figure
hist(cat(1,ETparams.saccadeInfo.duration),40)
xlabel('Saccade duration (s)'),ylabel('Number of saccades')
figure
hist(cat(1,ETparams.fixationInfo.duration),100)
xlabel('Fixation duration (s)'),ylabel('Number of fixations')