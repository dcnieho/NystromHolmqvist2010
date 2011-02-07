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
ETparams.screen.res.pix     = [1024 768];
ETparams.screen.size        = [0.38 0.30];
ETparams.screen.viewingDist = 0.67;

ETparams.samplingFreq = 1250;

ETparams.blink.velocityThreshold = 1000;            % if vel > 1000 degrees/s, it is noise or blinks
ETparams.blink.accThreshold = 100000;               % if acc > 100000 degrees/s^2, it is noise or blinks

ETparams.saccade.peakVelocityThreshold = 100;       % Initial value of the peak detection threshold.
ETparams.saccade.minDur = 10;                       % in milliseconds

ETparams.glissade.searchWindow = 40;                % window after saccade in which we search for glissades, in milliseconds
ETparams.glissade.maxDur = 80;                      % in milliseconds

ETparams.fixation.minDur = 40;                      % in milliseconds

% process params
[angleInPixelsH, angleInPixelsV] = degrees2pixels(1, ETparams.screen.viewingDist, ETparams.screen.res.pix, ETparams.screen.size);
ETparams.screen.res.deg     = [1024 768]./[angleInPixelsH angleInPixelsV];


%%-------------------------------------------------------------------------
%%% Begin detection
%%-------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Process eye-movement data for file (participant) and trial separately, i - files, j -
% trials
%--------------------------------------------------------------------------
data = cell(size(ETdata));
for i = 1:size(ETdata,1)
    for j = 1:size(ETdata,2)
        % Process data
        fprintf('Subj %d, Trial %d\n',i,j);
        data{i,j} = eventDetection(ETdata(i,j).X,ETdata(i,j).Y,ETparams);
        
        if 1
            % plot the trial (eye X, eye Y, velocity traces and scanpath,
            % as well as detected events (TODO: scanpath)
            fhndl = figure('Units','normalized','Position',[0 0 1 1]);  % make fullscreen figure
            plotDetection(data{i,j},ETparams.samplingFreq,ETparams.glissade.searchWindow,ETparams.screen.res.deg,sprintf('Subj %d, Trial %d',i,j));
            set(fhndl,'Visible','on');  % assert visibility to bring window to front again after keypress
            pause
            if ~ishghandle(fhndl)
                break;
            end
        end
    end
end

if 0
    save([cd,'\DetectionResults\DetectionResults.mat'],'data');
end

%%-------------------------------------------------------------------------
%%% Plot results
%%-------------------------------------------------------------------------

% Calculate basic parameters
catdata = [data{:}];
mean(cell2mat(arrayfun(@(x) x.duration,[catdata.glissade],'UniformOutput',false)))
mean(cell2mat(arrayfun(@(x) x.duration,[catdata.saccade ],'UniformOutput',false)))
mean(cell2mat(arrayfun(@(x) x.duration,[catdata.fixation],'UniformOutput',false)))

%Plot histograms
% figure
% hist(cat(1,ETparams.glissadeInfo.duration),100)
figure
hist(cell2mat(arrayfun(@(x) x.duration,[catdata.saccade ],'UniformOutput',false)),40)
xlabel('Saccade duration (s)'),ylabel('Number of saccades')
figure
hist(cell2mat(arrayfun(@(x) x.duration,[catdata.fixation],'UniformOutput',false)),100)
xlabel('Fixation duration (s)'),ylabel('Number of fixations')

plotDetection(data{2,4});