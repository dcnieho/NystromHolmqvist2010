function plotDetection(varargin)

% standard plot routine for monocular data
% See also plotWithMark, plot2D and subplot
% Two call syntaxes, one with all data in a struct, the other with all
% argument separately

if nargin<=2 && isstruct(varargin{1})
    %%% unpack the needed variables
    data = varargin{1};
    % time series
    xdata   = data.deg.X;
    ydata   = data.deg.Y;
    vel     = data.deg.vel;
    time    = ([1:length(xdata)]-1)/data.sampleRate * 1000;
    % markers
    sacon   = data.saccade.on;
    sacoff  = data.saccade.off;
    glisoff = data.glissade.off;
    glistyp = data.glissade.type;
    fixon   = data.fixation.on;
    fixoff  = data.fixation.off;
    % thresholds
    saccadePeakVelocityThreshold    = data.saccade.peakVelocityThreshold;
    saccadeOnsetVelocityTreshold    = data.saccade.onsetVelocityTreshold;
    saccadeOffsetVelocityTreshold   = data.saccade.offsetVelocityTreshold;
    glissadeSearchSamples           = ceil(data.glissadeSearchWindow/data.sampleRate * 1000);
    
    if nargin==1
        titel = '';
    else
        titel = varargin{2};
    end
end

% determine axis limits
mmxy = minMax2(xdata,ydata);
mmt  = [min(time) max(time)];

% plot X trace with fixation markers
ax = subplot('position',[0.10 0.84 0.80 0.12]);
plotWithMark(time,xdata,...                                             % data (y,x)
             'time (ms) - fixations','Horizontal (°)',titel,...         % y-axis label, x-axis label, axis title
             fixon, {'bo','MarkerFaceColor','blue','MarkerSize',4},...  % fixation on  markers
             fixoff,{'ro','MarkerFaceColor','red' ,'MarkerSize',4} ...  % fixation off markers
            );
axis([mmt(1) mmt(2) mmxy(1) mmxy(2)]);

% plot Y trace with fixation markers
ay = subplot('position',[0.10 0.68 0.80 0.12]);
plotWithMark(time,ydata,...                                             % data (y,x)
             'time (ms) - fixations','Vertical (°)','',...              % y-axis label, x-axis label, axis title
             fixon, {'bo','MarkerFaceColor','blue','MarkerSize',4},...  % fixation on  markers
             fixoff,{'ro','MarkerFaceColor','red' ,'MarkerSize',4} ...  % fixation off markers
            );
axis([mmt(1) mmt(2) mmxy(1) mmxy(2)]);

% plot velocity trace with saccade and glissade markers
qhighvelglissade = glistyp==2;      % determine glissade type: 1 is low velocity, 2 is high velocity
av = subplot('position',[0.10 0.52 0.80 0.12]);
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