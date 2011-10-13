function plotDetection(data,datatype,veltype,sampleRate,glissadeSearchWindow,rect,titel,pic)

% standard plot routine for monocular data
% See also plotWithMark, plot2D and subplot
% call syntax:
% plotDetection(data, datatype, sampleRate, glissadeSearchWindow, res, title)
% - data: see below for the fields it needs to contain, they're all
%         unpacked in the same place so its easy to see
% - datatype: 'pix' or 'deg'
% - veltype: type of velocity to plot. 'vel': 2D velocity,
%   'velX': X/azimuthal velocity, 'velY': Y/elevational velocity
% - sampleRate
% - glissadeSearchWindow: in milliseconds
% - rect: struct with extends of screen window, [upper left (x,y) lower
%   right (x,y)]. rect will be read from rect.(datatype)
% - title: can be empty
% - pic: optional. struct with two fields, imdata with the image and offset
%       to encode the offset between the top left of the screen and of the
%       picture.
%
% LEGEND of the plots:
% First two plots show the eye's azimuth and elevation in degree (Fick
% angles)
% Blue markers indicate the start of fixations, red the end
% Third plot shows the angular velocity of the eye, with blue markers
% indicating the start of saccades, red markers their end, cyan stars the
% end of high velocity glissades and green stars the end of low velocity
% glissades (start of glissades is end of preceding saccade)
% Finally, the last two plots show the subject's scanpath, either
% indicating the detected fixations or the raw eye position data. The eye
% position at the start of the trial (or the first fixation) is marked by a
% blue marker and the end of the trial (or last fixation) by a red marker.

error(nargchk(6,8,nargin,'struct'))

%%% unpack the needed variables
if nargin<7
    titel = '';
end

rect = rect.(datatype);

% prepare labels
switch veltype
    case 'vel'
        vlbltype = '2D';
    case 'velX'
        if strcmp(datatype,'deg')
            vlbltype = 'Azi';
        else
            vlbltype = 'X';
        end
    case 'velY'
        if strcmp(datatype,'deg')
            vlbltype = 'Ele';
        else
            vlbltype = 'Y';
        end
end
if strcmp(datatype,'deg')
    unit = '°';
    xlbl = ['Azimuth (' unit ')'];
    ylbl = ['Elevation (' unit ')'];
else
    unit = 'pix';
    xlbl = ['Horizontal (' unit ')'];
    ylbl = ['Vertical (' unit ')'];
end
vlbl = ['Velocity ' vlbltype ' (' unit '/s)'];
clbl = 'Xcorr  response';   % double space on purpose, reads easier for me

% time series
% position
if strcmp(datatype,'pix')
    xdata   = data.pix.X;
    ydata   = data.pix.Y;
elseif strcmp(datatype,'deg')
    xdata   = data.deg.Azi;
    ydata   = data.deg.Ele;
end

% time
time    = ([1:length(xdata)]-1)/sampleRate * 1000;

% velocity
if strcmp(datatype,'pix')
    vel     = data.pix.(veltype);
elseif strcmp(datatype,'deg')
    switch veltype
        case 'vel'
            vel     = data.deg.vel;
        case 'velX'
            vel     = data.deg.velAzi;
        case 'velY'
            vel     = data.deg.velEle;
    end
end

% markers
sacon   = data.saccade.on;
sacoff  = data.saccade.off;
if isfield(data,'blink')
    blinkMarks = {data.blink.on, {'mo','MarkerFaceColor','magenta','MarkerSize',4}, ... % blink on  markers
                  data.blink.off,{'mo','MarkerFaceColor','magenta','MarkerSize',4}};    % blink off markers
else
    blinkMarks = {};
end
if isfield(data,'glissade')
    qhighvelglissade = data.glissade.type==2;                                           % determine glissade type: 1 is low velocity, 2 is high velocity
    glismarks = {data.glissade.off(qhighvelglissade) ,{'c*'},...                        % high velocity glissade off markers
                 data.glissade.off(~qhighvelglissade),{'g*'}};                          % low  velocity glissade off markers
else
    glismarks = {};
end
if isfield(data,'fixation')
    qHaveFixations = true;
    xfixpos = data.fixation.(['meanX_' datatype]);
    yfixpos = data.fixation.(['meanY_' datatype]);
    fixmarks  = {data.fixation.on, {'bo','MarkerFaceColor','blue','MarkerSize',4},...   % fixation on  markers
                 data.fixation.off,{'ro','MarkerFaceColor','red' ,'MarkerSize',4}};     % fixation off markers
else
    qHaveFixations = false;
    fixmarks = {};
end
% thresholds
if isfield(data.saccade,'peakXCorrThreshold')
    % used saccade template
    qSaccadeTemplate = true;
    saccadePeakXCorrThreshold       = data.saccade.peakXCorrThreshold;
else
    % detected saccades based on velocity trace
    qSaccadeTemplate = false;
end

if isfield(data.saccade,'offsetXCorrThreshold')
    % refinement also run from xcorr responses
    qSaccadeTemplateRefinement      = true;
    saccadeOnsetXCorrThreshold      = data.saccade.onsetXCorrThreshold;
    saccadeOffsetXCorrThreshold     = data.saccade.offsetXCorrThreshold;
else
    % refinement run from velocity trace
    qSaccadeTemplateRefinement      = false;
    saccadePeakVelocityThreshold    = data.saccade.peakVelocityThreshold;
    saccadeOnsetVelocityThreshold   = data.saccade.onsetVelocityThreshold;
    saccadeOffsetVelocityThreshold  = data.saccade.offsetVelocityThreshold;
end
glissadeSearchSamples           = ceil(glissadeSearchWindow/sampleRate * 1000);


%%% determine time axis limits
mmt  = [min(time) max(time)];

%%% determine subplot positions
if qSaccadeTemplate
    xplotPos = [0.05 0.88 0.90 0.08];
    yplotPos = [0.05 0.76 0.90 0.08];
    vplotPos = [0.05 0.60 0.90 0.12];
    cplotPos = [0.05 0.44 0.90 0.12];
    fixplotPos = [0.05 0.06 0.43 0.34];
    rawplotPos = [0.52 0.06 0.43 0.34];
else
    xplotPos = [0.05 0.84 0.90 0.12];
    yplotPos = [0.05 0.68 0.90 0.12];
    vplotPos = [0.05 0.52 0.90 0.12];
    fixplotPos = [0.05 0.06 0.43 0.40];
    rawplotPos = [0.52 0.06 0.43 0.40];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% plot X trace with fixation markers
ax = subplot('position',xplotPos);
plotWithMark(time,xdata,...                                             % data (y,x)
             'time (ms) - fixations',xlbl,titel,...                     % x-axis label, y-axis label, axis title
             fixmarks{:} ...                                            % fixation markers (if any)
            );
axis([mmt(1) mmt(2) rect(1) rect(3)]);
axis ij


%%% plot Y trace with fixation markers
ay = subplot('position',yplotPos);
plotWithMark(time,ydata,...                                             % data (y,x)
             'time (ms) - fixations',ylbl,'',...                        % x-axis label, y-axis label, axis title
             fixmarks{:}, ...                                           % fixation markers (if any)
             blinkMarks{:} ...                                          % blink markers (if any)
            );
axis([mmt(1) mmt(2) rect(2) rect(4)]);
axis ij


%%% plot velocity trace with saccade and glissade markers
av = subplot('position',vplotPos);
% line at 0
plot([time(1) time(end)],[0 0],'b');
hold on;
plotWithMark(time,vel,...                                               % data (y,x)
             'time (ms) - saccades/glissades',vlbl,'',...               % x-axis label, y-axis label, axis title
             sacon, {'bo','MarkerFaceColor','blue','MarkerSize',4},...  % saccade on  markers
             sacoff,{'ro','MarkerFaceColor','red' ,'MarkerSize',4},...  % saccade off markers
             glismarks{:}, ...                                          % glissade markers (if any)
             blinkMarks{:} ...                                          % blink markers (if any)
            );
% add detection thresholds
if strcmp(datatype,'deg') && ~qSaccadeTemplateRefinement && strcmp(veltype,'vel')
    % dont plot if:
    % 1. if plotting pixels, as thresholds are in °/s
    % 2. if refinement was done with the saccade template responses as no
    %    velocity thresholds are then used; it would be misleading to plot
    %    them here
    % 3. if we're plotting a component velocity, as thresholds are for 2D
    %    velocity
    hold on;
    plot(mmt,[1 1]*saccadePeakVelocityThreshold,'r--')
    plot(mmt,[1 1]*saccadeOnsetVelocityThreshold,'r:')
    for p=1:length(sacoff)
        plot(time([sacoff(p) min(glissadeSearchSamples+sacoff(p),end)]),[1 1]*saccadeOffsetVelocityThreshold(p),'r-'); % the "end" is returns the length of time. Cool end works everywhere inside an index expression!
    end
    hold off;
end
axis([mmt(1) mmt(2) min(0,min(vel)) max(vel)]);
if ~strcmp(veltype,'vel')
    axis ij
end

%%% plot cross correlation output with saccade and glissade markers
if qSaccadeTemplate
    ac = subplot('position',cplotPos);
    % line at 0
    plot([time(1) time(end)],[0 0],'b');
    hold on;
    plotWithMark(time,data.deg.velXCorr,...                                 % data (y,x)
                 'time (ms) - saccades/glissades',clbl,'',...               % x-axis label, y-axis label, axis title
                 sacon, {'bo','MarkerFaceColor','blue','MarkerSize',4},...  % saccade on  markers
                 sacoff,{'ro','MarkerFaceColor','red' ,'MarkerSize',4},...  % saccade off markers
                 glismarks{:}, ...                                          % glissade markers (if any)
                 blinkMarks{:} ...                                          % blink markers (if any)
                );
    hold on;
    % add detection thresholds
    plot(mmt,[1 1]*saccadePeakXCorrThreshold,'r--')
    if qSaccadeTemplateRefinement
        plot(mmt,[1 1]*saccadeOnsetXCorrThreshold,'r:')
        for p=1:length(sacoff)
            plot(time([sacoff(p) min(glissadeSearchSamples+sacoff(p),end)]),[1 1]*saccadeOffsetXCorrThreshold(p),'r-'); % the "end" is returns the length of time. Cool end works everywhere inside an index expression!
        end
    end
    hold off;
    axis([mmt(1) mmt(2) 0 min(2.5,max(data.deg.velXCorr))]);    % xcorr values above 2.5 seem to only occur due to noise
else
    ac = [];
end

% link x-axis (time) of the three/four timeseries for easy viewing
linkaxes([ax ay av ac],'x');


%%% plot scanpath of raw data and of fixations
if qHaveFixations
    asf = subplot('position',fixplotPos);
    if nargin>=8 && strcmp(datatype,'pix') && ~isempty(pic)
        imagesc([0 size(pic.imdata,2)]+pic.offset(2),[0 size(pic.imdata,1)]+pic.offset(1),pic.imdata);
        hold on
    end
    plotWithMark(xfixpos,yfixpos,...                                                    % data (y,x)
                 xlbl,ylbl,'',...                                                       % x-axis label, y-axis label, axis title
                 [1:length(xfixpos)],{'go','MarkerFaceColor','g'   ,'MarkerSize',4},... % mark each fixation (that is marker on each datapoint we feed it
                 1,                  {'bo','MarkerFaceColor','blue','MarkerSize',4},... % make first fixation marker blue
                 length(xfixpos),    {'ro','MarkerFaceColor','red' ,'MarkerSize',4} ... % make last  fixation marker red
                );
    axis(rect([1 3 2 4]));
    axis ij
else
    asf = [];
end

asr = subplot('position',rawplotPos);
if nargin>=8 && strcmp(datatype,'pix') && ~isempty(pic)
    imagesc([0 size(pic.imdata,2)]+pic.offset(2),[0 size(pic.imdata,1)]+pic.offset(1),pic.imdata);
    hold on
end
plotWithMark(xdata,ydata,...                                                        % data (y,x)
             xlbl,ylbl,'',...                                                       % x-axis label, y-axis label, axis title
             1,                  {'bo','MarkerFaceColor','blue','MarkerSize',4},... % use blue marker for first datapoint
             length(xdata),      {'ro','MarkerFaceColor','red' ,'MarkerSize',4} ... % use red  marker for last  datapoint
            );
axis(rect([1 3 2 4]));
axis ij

% link view of the two scanpath plots for easy viewing
linkaxes([asr asf],'xy');


zoom on;