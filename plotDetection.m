function plotDetection(varargin)

% standard plot routine for monocular data
% See also plotWithMark, plot2D and subplot
% Two call syntaxes, one with most data in a struct, the other with all
% argument separately:
% plotDetection(data, sampleRate, glissadeSearchWindow, res, title)
% - data: see below for the fields it needs to contain, they're all
%         unpacked in the same place so its easy to see
% - sampleRate
% - glissadeSearchWindow: in milliseconds
% - rect: extends of screen window, [upper left (x,y) lower right (x,y)].
% - title: can be empty
% TODO: make way to select pixels or degrees for plotting (except velocity,
% which is always plotted in °/s

if nargin<=5 && isstruct(varargin{1})
    %%% unpack the needed variables
    if nargin==5
        titel = varargin{5};
    elseif nargin==4
        titel = '';
    else
        error('not enough input arguments');
    end
    
    data                 = varargin{1};
    sampleRate           = varargin{2};
    glissadeSearchWindow = varargin{3};
    rect                 = varargin{4};
    % time series
    xdata   = data.deg.X;
    ydata   = data.deg.Y;
    vel     = data.deg.vel;
    time    = ([1:length(xdata)]-1)/sampleRate * 1000;
    % markers
    sacon   = data.saccade.on;
    sacoff  = data.saccade.off;
    glisoff = data.glissade.off;
    glistyp = data.glissade.type;
    fixon   = data.fixation.on;
    fixoff  = data.fixation.off;
    xfixpos = data.fixation.meanX;
    yfixpos = data.fixation.meanY;
    % thresholds
    saccadePeakVelocityThreshold    = data.saccade.peakVelocityThreshold;
    saccadeOnsetVelocityTreshold    = data.saccade.onsetVelocityTreshold;
    saccadeOffsetVelocityTreshold   = data.saccade.offsetVelocityTreshold;
    glissadeSearchSamples           = ceil(glissadeSearchWindow/sampleRate * 1000);
else
    assert(nargin>34) %TODO
end


%%% determine time axis limits
mmt  = [min(time) max(time)];

%%% plot X trace with fixation markers
ax = subplot('position',[0.05 0.84 0.90 0.12]);
plotWithMark(time,xdata,...                                             % data (y,x)
             'time (ms) - fixations','Horizontal (°)',titel,...         % y-axis label, x-axis label, axis title
             fixon, {'bo','MarkerFaceColor','blue','MarkerSize',4},...  % fixation on  markers
             fixoff,{'ro','MarkerFaceColor','red' ,'MarkerSize',4} ...  % fixation off markers
            );
axis([mmt(1) mmt(2) rect(1) rect(3)]);

% plot Y trace with fixation markers
ay = subplot('position',[0.05 0.68 0.90 0.12]);
plotWithMark(time,ydata,...                                             % data (y,x)
             'time (ms) - fixations','Vertical (°)','',...              % y-axis label, x-axis label, axis title
             fixon, {'bo','MarkerFaceColor','blue','MarkerSize',4},...  % fixation on  markers
             fixoff,{'ro','MarkerFaceColor','red' ,'MarkerSize',4} ...  % fixation off markers
            );
axis([mmt(1) mmt(2) rect(2) rect(4)]);

% plot velocity trace with saccade and glissade markers
qhighvelglissade = glistyp==2;      % determine glissade type: 1 is low velocity, 2 is high velocity
av = subplot('position',[0.05 0.52 0.90 0.12]);
plotWithMark(time,vel,...                                               % data (y,x)
             'time (ms) - saccades/glissades','Velocity (°/s)','',...   % y-axis label, x-axis label, axis title
             sacon, {'bo','MarkerFaceColor','blue','MarkerSize',4},...  % saccade on  markers
             sacoff,{'ro','MarkerFaceColor','red' ,'MarkerSize',4},...  % saccade off markers
             glisoff(qhighvelglissade),{'c*'}                     ,...  % high velocity glissade off markers
             glisoff(~qhighvelglissade),{'g*'}                     ...  % low  velocity glissade off markers
            );
% add detection thresholds
hold on;
plot(mmt,[1 1]*saccadePeakVelocityThreshold,'r--')
plot(mmt,[1 1]*saccadeOnsetVelocityTreshold,'r:')
for p=1:length(sacoff)
    plot(time([sacoff(p) min(glissadeSearchSamples+sacoff(p),end)]),[1 1]*saccadeOffsetVelocityTreshold(p),'r-'); % the "end" is returns the length of time. Cool end works everywhere insde and index expression!
end
hold off;
axis([mmt(1) mmt(2) 0 max(vel)]);

% link x-axis (time) of the three timeseries for easy viewing
linkaxes([ax ay av],'x');


%%% plot scanpath of raw data and of fixations
asf = subplot('position',[0.05 0.06 0.43 0.40]);
plotWithMark(xfixpos,yfixpos,...                                                    % data (y,x)
             'Hor (°)','Ver (°)','',...                                             % y-axis label, x-axis label, axis title
             [1:length(xfixpos)],{'go','MarkerFaceColor','g'   ,'MarkerSize',4},... % mark each fixation (that is marker on each datapoint we feed it
             1,                  {'bo','MarkerFaceColor','blue','MarkerSize',4},... % make first fixation marker blue
             length(xfixpos),    {'ro','MarkerFaceColor','red' ,'MarkerSize',4} ... % make last  fixation marker red
            );
axis(rect([1 3 2 4]));
axis ij

asr = subplot('position',[0.52 0.06 0.43 0.40]);
plotWithMark(xdata,ydata,...                                                        % data (y,x)
             'Hor (°)','Ver (°)','',...                                             % y-axis label, x-axis label, axis title
             1,                  {'bo','MarkerFaceColor','blue','MarkerSize',4},... % use blue marker for first datapoint
             length(xdata),      {'ro','MarkerFaceColor','red' ,'MarkerSize',4} ... % use red  marker for last  datapoint
            );
axis(rect([1 3 2 4]));
axis ij

% link view of the two scanpath plots for easy viewing
linkaxes([asr asf],'xy');


zoom on;