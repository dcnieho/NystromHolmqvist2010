clear all, clear functions, close all; clc
%%-------------------------------------------------------------------------
%%% This code is an implementation of 
%%% Nyström, M. & Holmqvist, K. (in press), "An adaptive algorithm for
%%% fixation, saccade, and glissade detection in eye-tracking data".
%%% Behavior Research Methods
%%-------------------------------------------------------------------------

%%-------------------------------------------------------------------------
%%% Init parameters
%%-------------------------------------------------------------------------
load('1250Hz_3_Participants.mat');

% load parameters
ETparams = getParameters;

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
            % as well as detected events)
            if ~ishghandle(fhndl)
                fhndl = figure('Units','normalized','Position',[0 0 1 1]);  % make fullscreen figure
            end
            plotDetection(data{i,j},'deg',ETparams.samplingFreq,ETparams.glissade.searchWindow,ETparams.screen.rect,sprintf('Subj %d, Trial %d',i,j));
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