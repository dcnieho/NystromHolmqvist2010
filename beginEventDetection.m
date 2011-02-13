clear all, clear functions, close all; clc

%%-------------------------------------------------------------------------
%%% Init parameters
%%-------------------------------------------------------------------------
load('1250Hz_3_Participants.mat');

% user settings
ETparams.screen.resolution              = [1024 768];
ETparams.screen.size                    = [0.38 0.30];
ETparams.screen.viewingDist             = 0.67;
ETparams.screen.subjectStraightAhead    = [512 200];    % Specify the screen coordinate that is straight ahead of the subject. Just specify the middle of the screen unless its important to you to get this very accurate!

% flip the Y coordinate of the data? All the routines assume the origin of
% the screen (0,0) is at the top left corner. You'll have to flip if the
% your data's origin is the lower left corner. Do a flip X if your origin
% is on the right side of the screen (sic).
ETparams.data.qFlipY                    = false;
ETparams.data.qFlipX                    = false;
% Do a precise calculation of angular eye velocity and acceleration? If
% not, we compute derivatives of eye azimuth and elevation analytically
% from the parameters of a fitted polynomial and then apply Pythagoras'
% theorem to compute eye velocity/acceleration. This is crude and should
% not be used if you're interested in the eye velocity, but its sufficient
% if you simply want to detect saccades in periods of fixation and/or
% smooth pursuit and are not interested in accurate measures of eye
% velocity/acceleration.
ETparams.data.qPreciseCalcDeriv         = false;

ETparams.samplingFreq                   = 1250;

ETparams.blink.velocityThreshold        = 1000;         % if vel > 1000 �/s, it is noise or blinks
ETparams.blink.accThreshold             = 100000;       % if acc > 100000 �/s�, it is noise or blinks

ETparams.saccade.peakVelocityThreshold  = 100;          % Initial value of the peak detection threshold, �/s
ETparams.saccade.minDur                 = 10;           % in milliseconds

ETparams.glissade.searchWindow          = 40;           % window after saccade in which we search for glissades, in milliseconds
ETparams.glissade.maxDur                = 80;           % in milliseconds

ETparams.fixation.minDur                = 40;           % in milliseconds
% How to deal with NaNs during possible fixation periods:
% 1: do not allow NaN during fixations
% 2: ignore NaNs and calculate mean fixation position based on other data
%    (not recommended in almost any situation, if you don't like 1,
%    consider option 3)
% 3: split fixation into multiple, providing each is at least minDur long
%    (e.g. one 250 ms fixation with some data missing in the middle might
%    be split up into a 100 ms and a 120 ms fixation)
ETparams.fixations.treatNaN             = 3;

% process params
ETparams = prepareParameters(ETparams);


%%-------------------------------------------------------------------------
%%% Begin detection
%%-------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Process eye-movement data, per participant (i), per trial (j)
%--------------------------------------------------------------------------
data = cell(size(ETdata));
fhndl = -1;
for i = 2%1:size(ETdata,1)
    for j = 6%1:size(ETdata,2)
        % Process data
        fprintf('Subj %d, Trial %d\n',i,j);
        data{i,j} = eventDetection(ETdata(i,j).X,ETdata(i,j).Y,ETparams);
        
        if 1
            % plot the trial (eye X, eye Y, velocity traces and scanpath,
            % as well as detected events
            if ~ishghandle(fhndl)
                fhndl = figure('Units','normalized','Position',[0 0 1 1]);  % make fullscreen figure
            end
            plotDetection(data{i,j},ETparams.samplingFreq,ETparams.glissade.searchWindow,ETparams.screen.rect.deg,sprintf('Subj %d, Trial %d',i,j));
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

figure('Units','normalized','Position',[0 0 1 1]);  % make fullscreen figure
plotDetection(data{2,4},ETparams.samplingFreq,ETparams.glissade.searchWindow,ETparams.screen.rect.deg,sprintf('Subj %d, Trial %d',2,4));