function plotDetectionFilt(data,datatype,sampleRate,glissadeSearchWindow,rect,titel,pic)

% standard plot routine for monocular data
% See also plotWithMark, plot2D and subplot
% call syntax:
% plotDetection(data, datatype, sampleRate, glissadeSearchWindow, res, title)
% - data: see below for the fields it needs to contain, they're all
%         unpacked in the same place so its easy to see
% - datatype: 'pix' or 'deg'
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

error(nargchk(5,7,nargin,'struct'))

%%% unpack the needed variables
if nargin<6
    titel = '';
end

rect = rect.(datatype);
if strcmp(datatype,'deg')
    unit = '°';
    xlbl = ['Azimuth (' unit ')'];
    ylbl = ['Elevation (' unit ')'];
    vlbl = ['Velocity (' unit '/s)'];
else
    unit = 'pix';
    xlbl = ['Horizontal (' unit ')'];
    ylbl = ['Vertical (' unit ')'];
    vlbl = ['Velocity (' unit '/s)'];
end

% time series
xdata   = data.(datatype).X;
ydata   = data.(datatype).Y;
vel     = data.(datatype).vel;

assert(isfield(data.(datatype),'velfilt'),'Saccades have not been cut from this data')
velcut  = data.(datatype).velfilt;
if isfield(data.(datatype),'Xfilt')
    xdatcut = data.(datatype).Xfilt;
    ydatcut = data.(datatype).Yfilt;
    qReconstructPos = true;
else
    qReconstructPos = false;
end
time    = ([1:length(xdata)]-1)/sampleRate * 1000;
% markers
sacon   = data.saccade.on;
sacoff  = data.saccade.off;
if isfield(data,'glissade')
    qHaveGlissades = true;
    glisoff = data.glissade.off;
    glistyp = data.glissade.type;
else
    qHaveGlissades = false;
end
if isfield(data,'fixation')
    qHaveFixations = true;
    fixon   = data.fixation.on;
    fixoff  = data.fixation.off;
    xfixpos = data.fixation.(['meanX_' datatype]);
    yfixpos = data.fixation.(['meanY_' datatype]);
else
    qHaveFixations = false;
end
% thresholds
saccadePeakVelocityThreshold    = data.saccade.peakVelocityThreshold;
saccadeOnsetVelocityTreshold    = data.saccade.onsetVelocityTreshold;
saccadeOffsetVelocityTreshold   = data.saccade.offsetVelocityTreshold;
glissadeSearchSamples           = ceil(glissadeSearchWindow/sampleRate * 1000);


%%% determine time axis limits
mmt  = [min(time) max(time)];

%%% plot X trace with fixation markers
ax = subplot('position',[0.05 0.84 0.90 0.12]);
if qHaveFixations
    fixmarks = {fixon, {'bo','MarkerFaceColor','blue','MarkerSize',4},...  % fixation on  markers
                fixoff,{'ro','MarkerFaceColor','red' ,'MarkerSize',4}};    % fixation off markers
else
    fixmarks = {};
end
hold on;
if qReconstructPos
    plot(time,xdata,'g');
    plotWithMark(time,xdatcut,...                                           % data (y,x)
                 'time (ms) - fixations',xlbl,titel,...                     % x-axis label, y-axis label, axis title
                 fixmarks{:} ...                                            % fixation markers (if any)
                );
else
    plotWithMark(time,xdata,...                                             % data (y,x)
                 'time (ms) - fixations',xlbl,titel,...                     % x-axis label, y-axis label, axis title
                 fixmarks{:} ...                                            % fixation markers (if any)
                );
end
axis([mmt(1) mmt(2) rect(1) rect(3)]);
axis ij


%%% plot Y trace with fixation markers
ay = subplot('position',[0.05 0.68 0.90 0.12]);
hold on;
if qReconstructPos
    plot(time,ydata,'g');
    plotWithMark(time,ydatcut,...                                           % data (y,x)
                 'time (ms) - fixations',ylbl,'',...                        % x-axis label, y-axis label, axis title
                 fixmarks{:} ...                                            % fixation markers (if any)
                );
else
    plotWithMark(time,ydata,...                                           % data (y,x)
                 'time (ms) - fixations',ylbl,'',...                        % x-axis label, y-axis label, axis title
                 fixmarks{:} ...                                            % fixation markers (if any)
                );
end
axis([mmt(1) mmt(2) rect(2) rect(4)]);
axis ij


%%% plot velocity trace with saccade and glissade markers
av = subplot('position',[0.05 0.52 0.90 0.12]);
if qHaveGlissades
    qhighvelglissade = glistyp==2;                                      % determine glissade type: 1 is low velocity, 2 is high velocity
    glismarks = {glisoff(qhighvelglissade) ,{'c*'},...                  % high velocity glissade off markers
                 glisoff(~qhighvelglissade),{'g*'}};                    % low  velocity glissade off markers
else
    glismarks = {};
end
% at line at 0
plot([time(1) time(end)],[0 0],'b');
hold on;
plot(time,vel,'g');
plotWithMark(time,velcut,...                                            % data (y,x)
             'time (ms) - saccades/glissades',vlbl,'',...               % x-axis label, y-axis label, axis title
             sacon, {'bo','MarkerFaceColor','blue','MarkerSize',4},...  % saccade on  markers
             sacoff,{'ro','MarkerFaceColor','red' ,'MarkerSize',4},...  % saccade off markers
             glismarks{:} ...                                           % glissade markers (if any)
            );
hold on;
% add detection thresholds
if strcmp(datatype,'deg')
    % don't add if plotting pixels, it doesn't make any sense then
    hold on;
    plot(mmt,[1 1]*saccadePeakVelocityThreshold,'r--')
    plot(mmt,[1 1]*saccadeOnsetVelocityTreshold,'r:')
    for p=1:length(sacoff)
        plot(time([sacoff(p) min(glissadeSearchSamples+sacoff(p),end)]),[1 1]*saccadeOffsetVelocityTreshold(p),'r-'); % the "end" is returns the length of time. Cool end works everywhere inside an index expression!
    end
    hold off;
end
axis([mmt(1) mmt(2) 0 max(vel)]);

% link x-axis (time) of the three timeseries for easy viewing
linkaxes([ax ay av],'x');


%%% plot scanpath of raw data and of data with saccades cut out
asf = subplot('position',[0.05 0.06 0.43 0.40]);
if qReconstructPos
    if nargin>=7 && strcmp(datatype,'pix') && ~isempty(pic)
        imagesc([0 size(pic.imdata,2)]+pic.offset(2),[0 size(pic.imdata,1)]+pic.offset(1),pic.imdata);
        hold on
    end
    plotWithMark(xdatcut,ydatcut,...                                                    % data (y,x)
                 xlbl,ylbl,'',...                                                       % x-axis label, y-axis label, axis title
                 1,                  {'bo','MarkerFaceColor','blue','MarkerSize',4},... % use blue marker for first datapoint
                 length(xdata),      {'ro','MarkerFaceColor','red' ,'MarkerSize',4} ... % use red  marker for last  datapoint
                );
    axis(rect([1 3 2 4]));
    axis ij
end

asr = subplot('position',[0.52 0.06 0.43 0.40]);
if nargin>=7 && strcmp(datatype,'pix') && ~isempty(pic)
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
