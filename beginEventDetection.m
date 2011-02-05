%--------------------------------------------------------------------------
% For usage, see the README-file.
%--------------------------------------------------------------------------
clear all, close all, clc
global ETparams

%--------------------------------------------------------------------------
% Init parameters
%--------------------------------------------------------------------------
load('1250Hz_3_Participants.mat');

ETparams.data = ETdata;
ETparams.screenSz = [1024 768];
ETparams.screenDim = [0.38 0.30];
ETparams.viewingDist = 0.67;
ETparams.samplingFreq = 1250;
ETparams.blinkVelocityThreshold = 1000;             % if vel > 1000 degrees/s, it is noise or blinks
ETparams.blinkAccThreshold = 100000;               % if acc > 100000 degrees/s^2, it is noise or blinks
ETparams.peakDetectionThreshold = 100;              % Initial value of the peak detection threshold. 

ETparams.minFixDur = 0.040; % in seconds
ETparams.minSaccadeDur = 0.010; % in seconds

%--------------------------------------------------------------------------
% Begin detection
%--------------------------------------------------------------------------

% Process data
eventDetection

save([cd,'\DetectionResults\DetectionResults.mat'],'ETparams');

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------

% Calculate basic parameters
mean(cat(1,ETparams.data.avgNoise))
mean(cat(1,ETparams.data.stdNoise))
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


plotResultsVel(ETparams,2,4)