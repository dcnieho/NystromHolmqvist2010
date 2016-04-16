function plotDetectionFilt(data,datatype,veltype,sampleRate,glissadeSearchWindow,rect,titel,pic)

% standard plot routine for monocular data
% See also plotWithMark, plot2D and axes
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

narginchk(6,8)

assert(isfield(data.(datatype),'vel'),'data for %s not available',datatype);

%%% unpack the needed variables
if nargin<7
    titel = '';
else
    titel = texlabel(titel,'literal');
end

rect = rect.(datatype);

% prepare labels
missing = data.deg.missing;
if strcmp(datatype,'deg')
    unit = '°';
    xlbl = ['Azimuth (' unit ')'];
    ylbl = ['Elevation (' unit ')'];
else
    unit = 'pix';
    xlbl = ['Horizontal (' unit ')'];
    ylbl = ['Vertical (' unit ')'];
end
plbl = 'pupil size';
pvlbl= 'abs({\delta} pupil size)';

if strcmp(datatype,'deg')
    vxlbl = 'Azi';
    vylbl = 'Ele';
else
    vxlbl = 'X';
    vylbl = 'Y';
end
vlbl = {['Velocity ' '2D'  ' (' unit '/s)'],['Velocity ' vxlbl ' (' unit '/s)'], ['Velocity ' vylbl ' (' unit '/s)']};
vidx = find(ismember({'vel','velX','velY'},veltype));
clbl = 'Xcorr  response';   % double space on purpose, reads easier for me

% time series
assert(isfield(data.(datatype),'velFilt'),'Saccades have not been cut from this data')
% position
if strcmp(datatype,'pix')
    xdata   = data.pix.X;
    ydata   = data.pix.Y;
elseif strcmp(datatype,'deg')
    xdata   = data.deg.Azi;
    ydata   = data.deg.Ele;
end

if strcmp(datatype,'pix') && isfield(data.pix,'Xfilt')
    xdatcut = data.pix.XFilt;
    ydatcut = data.pix.YFilt;
    qReconstructPos = true;
elseif strcmp(datatype,'deg') && isfield(data.deg,'AziFilt')
    xdatcut = data.deg.AziFilt;
    ydatcut = data.deg.EleFilt;
    qReconstructPos = true;
else
    qReconstructPos = false;
end

% see if also have saccadic traces
if strcmp(datatype,'pix') && isfield(data.pix,'XSac')
    xSac    = data.pix.XSac;
    ySac    = data.pix.YSac;
    qHaveSaccadic = true;
elseif strcmp(datatype,'deg') && isfield(data.deg,'AziSac')
    xSac    = data.deg.AziSac;
    ySac    = data.deg.EleSac;
    qHaveSaccadic = true;
else
    qHaveSaccadic = false;
end

% time
time    = ([1:length(xdata)]-1)/sampleRate * 1000;

% velocity
if strcmp(datatype,'pix')
    vel     = {data.pix.vel,data.pix.velX,data.pix.velY};
    velcut  = {data.pix.velFilt,data.pix.velXFilt,data.pix.velYFilt};
    if qHaveSaccadic
        velSac  = {data.pix.velSac,data.pix.velXSac,data.pix.velYSac};
    else
        velSac  = {[],[],[]};
    end
elseif strcmp(datatype,'deg')
    vel     = {data.deg.vel,data.deg.velAzi,data.deg.velEle};
    velcut  = {data.deg.velFilt,data.deg.velAziFilt,data.deg.velEleFilt};
    if qHaveSaccadic
        velSac  = {data.deg.velSac,data.deg.velAziSac,data.deg.velEleSac};
    else
        velSac  = {[],[],[]};
    end
end

% markers
if ~isempty(missing.on)
    missFlag = Interleave(arrayfun(@(on,off) on:off,missing.on,missing.off,'uni',false),repmat({{'-r'}},1,length(missing.on)));
else
    missFlag = {};
end

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
    glisMarks = {data.glissade.off(qhighvelglissade) ,{'c*'},...                        % high velocity glissade off markers
                 data.glissade.off(~qhighvelglissade),{'g*'}};                          % low  velocity glissade off markers
else
    glisMarks = {};
end
if isfield(data,'fixation')
    fixMarks  = {data.fixation.on, {'bo','MarkerFaceColor','blue','MarkerSize',4},...   % fixation on  markers
                 data.fixation.off,{'ro','MarkerFaceColor','red' ,'MarkerSize',4}};     % fixation off markers
else
    fixMarks  = {};
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
glissadeSearchSamples   = ceil(glissadeSearchWindow./1000 * sampleRate);


%%% determine time axis limits
mmt  = [min(time) max(time)];

%%% determine axes positions
if qSaccadeTemplate
    if isfield(data,'pupil') && ~isempty(data.pupil.size)
        xplotPos = [0.05 0.88 0.90 0.08];
        yplotPos = [0.05 0.76 0.90 0.08];
        pplotPos = [0.05 0.64 0.90 0.08];
        vplotPos = [0.05 0.50 0.90 0.10];
        cplotPos = [0.05 0.36 0.90 0.10];
        cutplotPos = [0.05 0.04 0.43 0.28];
        rawplotPos = [0.52 0.04 0.43 0.28];
    else
        xplotPos = [0.05 0.88 0.90 0.08];
        yplotPos = [0.05 0.76 0.90 0.08];
        vplotPos = [0.05 0.60 0.90 0.12];
        cplotPos = [0.05 0.44 0.90 0.12];
        cutplotPos = [0.05 0.06 0.43 0.34];
        rawplotPos = [0.52 0.06 0.43 0.34];
    end
else
    if isfield(data,'pupil') && ~isempty(data.pupil.size)
        xplotPos = [0.05 0.88 0.90 0.08];
        yplotPos = [0.05 0.76 0.90 0.08];
        pplotPos = [0.05 0.60 0.90 0.12];
        vplotPos = [0.05 0.44 0.90 0.12];
        cutplotPos = [0.05 0.06 0.43 0.34];
        rawplotPos = [0.52 0.06 0.43 0.34];
    else
        xplotPos = [0.05 0.84 0.90 0.12];
        yplotPos = [0.05 0.68 0.90 0.12];
        vplotPos = [0.05 0.52 0.90 0.12];
        cutplotPos = [0.05 0.06 0.43 0.40];
        rawplotPos = [0.52 0.06 0.43 0.40];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% plot X trace with fixation markers
ax = axes('position',xplotPos);
hold on;
if qReconstructPos
    plot(time,xdata,'k');
    pdat = xdatcut;
    style = {'g-'};
else
    pdat = xdata;
    style = {'k-'};
end
if qHaveSaccadic
    plot(time,xSac,'c');
end
plotWithMark(time,pdat,style,...                                        % data (y,x), style
             'time (ms) - fixations',xlbl,titel,...                     % x-axis label, y-axis label, axis title
             missFlag{:}, ...                                           % color part of trace that is missing
             fixMarks{:}, ...                                           % fixation markers (if any)
             blinkMarks{:} ...                                          % blink markers (if any)
            );
axis([mmt(1) mmt(2) rect(1) rect(3)]);
axis ij


%%% plot Y trace with fixation markers
ay = axes('position',yplotPos);
hold on;
if qReconstructPos
    plot(time,ydata,'k');
    pdat = ydatcut;
    style = {'g-'};
else
    pdat = ydata;
    style = {'k-'};
end
if qHaveSaccadic
    plot(time,ySac,'c');
end
plotWithMark(time,pdat,style,...                                        % data (y,x), style
             'time (ms) - fixations',ylbl,'',...                        % x-axis label, y-axis label, axis title
             missFlag{:}, ...                                           % color part of trace that is missing
             fixMarks{:}, ...                                           % fixation markers (if any)
             blinkMarks{:} ...                                          % blink markers (if any)
            );
axis([mmt(1) mmt(2) rect(2) rect(4)]);
axis ij


%%% plot pupil size trace with blink markers
if isfield(data,'pupil') && ~isempty(data.pupil.size)
    % size
    ap = axes('position',pplotPos);
    plotWithMark(time,data.pupil.size,{'k-'},...                        % data (y,x), style
                 'time (ms) - blinks',plbl,'',...                       % x-axis label, y-axis label, axis title
                 missFlag{:}, ...                                       % color part of trace that is missing
                 blinkMarks{:} ...                                      % blink markers (if any)
                );
    axis([mmt(1) mmt(2) min(data.pupil.size) max(data.pupil.size)]);
    % change of size
    pvdat = abs(data.pupil.dsize);
    apv= axes('position',pplotPos);
    % line at 0
    plot([time(1) time(end)],[0 0],'b');
    hold on;
    plotWithMark(time,pvdat,{'k-'},...                                  % data (y,x), style
                 'time (ms) - blinks',pvlbl,'',...                      % x-axis label, y-axis label, axis title
                 missFlag{:}, ...                                       % color part of trace that is missing
                 blinkMarks{:} ...                                      % blink markers (if any)
                );
    axis([mmt(1) mmt(2) 0 max(pvdat)]);
    if isfield(data.blink,'peakDSizeThreshold')
        % plot pupil size change thresholds for blink detection
        hold on;
        plot(mmt, [1 1]*data.blink.peakDSizeThreshold,'r--')
        plot(mmt, [1 1]*data.blink.onsetDSizeThreshold,'r:')
        for p=1:length(data.blink.off)
            plot(time([data.blink.off(p) min(glissadeSearchSamples+data.blink.off(p),end)]), [1 1]*data.blink.offsetDSizeThreshold(p),'r-'); % the "end" is returns the length of time. Cool end works everywhere inside an index expression!
        end
        hold off;
    end
    % at start size, not dsize, is visible
    set([apv; allchild(apv)],'visible','off');
    % need to set visible axis to current and topmost axis for zooming to
    % operate on rigth axis. (If we don't do this, x zoom still works as
    % all axes are linked on x, by looking at y range reveals last drawn in
    % a given position is always target of zoom, even if it isn't
    % visible...)
    axes(ap);
    % toggle button
    uicontrol(...
    'Style','togglebutton',...
    'String','d',...
    'FontName','symbol',...
    'Units','Normalized',...
    'Position',[sum(pplotPos([1 3]))+.01 sum(pplotPos([2 4]))-.04 .02 .03],...
    'Callback',@Pupil_Callback);
else
    ap  = [];
    apv = [];
end


%%% plot velocity trace with saccade and glissade markers
av2 = axes('position',vplotPos);
plotVel(time,vel{1},velcut{1},velSac{1},vlbl{1},'vel',datatype,...
    missFlag,sacon,sacoff,glisMarks,blinkMarks,mmt,...
    qSaccadeTemplateRefinement,saccadePeakVelocityThreshold,saccadeOnsetVelocityThreshold,glissadeSearchSamples,saccadeOffsetVelocityThreshold);
avx = axes('position',vplotPos);
plotVel(time,vel{2},velcut{2},velSac{2},vlbl{2},'velX',datatype,...
    missFlag,sacon,sacoff,glisMarks,blinkMarks,mmt,...
    qSaccadeTemplateRefinement,saccadePeakVelocityThreshold,saccadeOnsetVelocityThreshold,glissadeSearchSamples,saccadeOffsetVelocityThreshold);
avy = axes('position',vplotPos);
plotVel(time,vel{3},velcut{3},velSac{3},vlbl{3},'velY',datatype,...
    missFlag,sacon,sacoff,glisMarks,blinkMarks,mmt,...
    qSaccadeTemplateRefinement,saccadePeakVelocityThreshold,saccadeOnsetVelocityThreshold,glissadeSearchSamples,saccadeOffsetVelocityThreshold);
vaxs = [av2 avx avy];
% show desired vel at start
toHide = [1:3]; toHide(toHide==vidx) = [];
for p=toHide
    set([vaxs(p); allchild(vaxs(p))],'visible','off');
end
axes(vaxs(vidx));   % set visible axis to current and topmost axis
% toggle button
strs = {'v2','vx','vy'};
strs2= strs(toHide);
vt1 = uicontrol(...
    'Style','pushbutton',...
    'String',strs2{1},...
    'Units','Normalized',...
    'Position',[sum(vplotPos([1 3]))+.01 sum(vplotPos([2 4]))-.04 .02 .03],...
    'Callback',@Velocity_Callback);
vt2 = uicontrol(...
    'Style','pushbutton',...
    'String',strs2{2},...
    'Units','Normalized',...
    'Position',[sum(vplotPos([1 3]))+.01 sum(vplotPos([2 4]))-.08 .02 .03],...
    'Callback',@Velocity_Callback);

%%% plot cross correlation output with saccade and glissade markers
if qSaccadeTemplate
    ac = axes('position',cplotPos);
    % line at 0
    plot([time(1) time(end)],[0 0],'b');
    hold on;
    plotWithMark(time,data.deg.velXCorr,{'k-'},...                          % data (y,x), style
                 'time (ms) - saccades/glissades',clbl,'',...               % x-axis label, y-axis label, axis title
                 sacon, {'bo','MarkerFaceColor','blue','MarkerSize',4},...  % saccade on  markers
                 sacoff,{'ro','MarkerFaceColor','red' ,'MarkerSize',4},...  % saccade off markers
                 glisMarks{:}, ...                                          % glissade markers (if any)
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
linkaxes([ax ay ap apv av2 avx avy ac],'x');


%%% plot scanpath of raw data and of data with saccades cut out
if qReconstructPos
    asf = axes('position',cutplotPos);
    if qHaveSaccadic
        plot(xSac,ySac,'c');
        hold on;
    end
    if nargin>=8 && strcmp(datatype,'pix') && ~isempty(pic)
        imagesc([0 size(pic.imdata,2)]+pic.offset(2),[0 size(pic.imdata,1)]+pic.offset(1),pic.imdata);
        hold on
    end
    plotWithMark(xdatcut,ydatcut,{'g-'},...                                             % data (y,x), style
                 xlbl,ylbl,'',...                                                       % x-axis label, y-axis label, axis title
                 1,                  {'bo','MarkerFaceColor','blue','MarkerSize',4},... % use blue marker for first datapoint
                 length(xdata),      {'ro','MarkerFaceColor','red' ,'MarkerSize',4} ... % use red  marker for last  datapoint
                );
    axis(rect([1 3 2 4]));
    axis ij
else
    asf = [];
end

asr = axes('position',rawplotPos);
if nargin>=8 && strcmp(datatype,'pix') && ~isempty(pic)
    imagesc([0 size(pic.imdata,2)]+pic.offset(2),[0 size(pic.imdata,1)]+pic.offset(1),pic.imdata);
    hold on
end
plotWithMark(xdata,ydata,{'k-'},...                                                 %  data (y,x), style
             xlbl,ylbl,'',...                                                       % x-axis label, y-axis label, axis title
             1,                  {'bo','MarkerFaceColor','blue','MarkerSize',4},... % use blue marker for first datapoint
             length(xdata),      {'ro','MarkerFaceColor','red' ,'MarkerSize',4} ... % use red  marker for last  datapoint
            );
axis(rect([1 3 2 4]));
axis ij

% link view of the two scanpath plots for easy viewing
linkaxes([asr asf],'xy');

% make sure we don't lose the standard toolbar
set(gcf,'Toolbar','figure');
zoom on;




    function Pupil_Callback(~,~,~)
        if strcmp(get(ap,'visible'),'on')
            set([ap;   allchild(ap)],'visible','off');
            set([apv; allchild(apv)],'visible','on');
            axes(apv);  % set visible axis to current and topmost axis
        else
            set([apv; allchild(apv)],'visible','off');
            set([ap;   allchild(ap)],'visible','on');
            axes(ap);   % set visible axis to current and topmost axis
        end
    end

    function Velocity_Callback(obj,~,~)
        % find which pressed, and therefore which to show and hide
        but = get(obj,'String');
        qShow  = strcmp(but,strs);
        toHidev = find(~qShow);
        toShow = find( qShow);
        % do show and hide
        for q=toHidev
            set([vaxs(q); allchild(vaxs(q))],'visible','off');
        end
        set([vaxs(toShow); allchild(vaxs(toShow))],'visible','on');
        axes(vaxs(toShow));     % set visible axis to current and topmost axis
        % adjust labels on buttons
        strs2v  = strs(toHidev);
        set(vt1,'String',strs2v{1});
        set(vt2,'String',strs2v{2});
    end

end

function plotVel(time,vel,velcut,velSac,vlbl,veltype,datatype,...
    missFlag,sacon,sacoff,glisMarks,blinkMarks,mmt,...
    qSaccadeTemplateRefinement,saccadePeakVelocityThreshold,saccadeOnsetVelocityThreshold,glissadeSearchSamples,saccadeOffsetVelocityThreshold)
% line at 0
plot([time(1) time(end)],[0 0],'b');
hold on;
plot(time,vel,'k');
if ~isempty(velSac)
    plot(time,velSac,'c');
end
plotWithMark(time,velcut,{'g-'},...                                     % data (y,x), style
             'time (ms) - saccades/glissades',vlbl,'',...               % x-axis label, y-axis label, axis title
             missFlag{:}, ...                                           % color part of trace that is missing
             sacon, {'bo','MarkerFaceColor','blue','MarkerSize',4},...  % saccade on  markers
             sacoff,{'ro','MarkerFaceColor','red' ,'MarkerSize',4},...  % saccade off markers
             glisMarks{:}, ...                                          % glissade markers (if any)
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
end