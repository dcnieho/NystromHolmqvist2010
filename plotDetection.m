function plotDetection(data,datatype,veltype,sampleRate,glissadeSearchWindow,rect,varargin)

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
% the rest is optional key-value parameters:
% - 'title': plot title
% - 'pic': struct with two fields, imdata with the image and offset to
%   encode the offset between the top left of the screen and of the
%   picture.
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

narginchk(6,inf)

assert(isfield(data.(datatype),'vel'),'data for %s not available',datatype);

%%% unpack the needed variables
% key-val parameters
titel = '';
pic = [];
if nargin>=7
    nKeyValInp = nargin-6;
    assert(mod(nKeyValInp,2)==0,'key-value arguments must come in pairs')
    expectVal = false;                          % start expecting an option name, not a value
    p = 1;
    while p <= length(varargin)
        if ~expectVal   % the current value should be a setting
            assert(ischar(varargin{p}),'option name must be a string')
            expectVal = true;
        else    % we just read a setting name, now look for a value for that setting
            switch varargin{p-1}
                case 'title'
                    titel = texlabel(varargin{p},'literal');
                case 'pic'
                    pic = varargin{p};
                case 'highlight'
                    if ~any(isnan(varargin{p}))
                        highlightTime = varargin{p};
                    end
                otherwise
                    error('do not understand input %s',varargin{p-1})
            end
            expectVal = false;
        end
        p=p+1;
    end
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
vlbl = {['Velocity 2D'      ' (' unit '/s)'],['Velocity '     vxlbl ' (' unit '/s)'], ['Velocity '     vylbl ' (' unit '/s)']};
albl = {['Acceleration 2D'  ' (' unit '/s)'],['Acceleration ' vxlbl ' (' unit '/s)'], ['Acceleration ' vylbl ' (' unit '/s)']};
vidx = find(ismember({'vel','velX','velY'},veltype));
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
    vel     = {data.pix.vel,data.pix.velX,data.pix.velY};
elseif strcmp(datatype,'deg')
    vel     = {data.deg.vel,data.deg.velAzi,data.deg.velEle};
end
% acceleration
if isfield(data.(datatype),'acc')
    qHaveAcceleration = true;
    if strcmp(datatype,'pix')
        acc     = {data.pix.acc,data.pix.accX,data.pix.accY};
    elseif strcmp(datatype,'deg')
        acc     = {data.deg.acc,data.deg.accAzi,data.deg.accEle};
    end
else
    qHaveAcceleration = false;
end

% markers
% for missing flags, also include blinks. we'd want to color original or
% interpolated data during a blink as well
if isfield(data,'blink')
    qMissOrBlink =                bounds2bool(missing.on   ,missing.off   ,length(vel{1}));
    qMissOrBlink = qMissOrBlink | bounds2bool(data.blink.on,data.blink.off,length(vel{1}));
    [missing.on,missing.off] = bool2bounds(qMissOrBlink);
end
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
        fixplotPos = [0.05 0.04 0.43 0.28];
        rawplotPos = [0.52 0.04 0.43 0.28];
    else
        xplotPos = [0.05 0.88 0.90 0.08];
        yplotPos = [0.05 0.76 0.90 0.08];
        vplotPos = [0.05 0.60 0.90 0.12];
        cplotPos = [0.05 0.44 0.90 0.12];
        fixplotPos = [0.05 0.06 0.43 0.34];
        rawplotPos = [0.52 0.06 0.43 0.34];
    end
elseif qHaveAcceleration
    if isfield(data,'pupil') && ~isempty(data.pupil.size)
        xplotPos = [0.05 0.88 0.90 0.08];
        yplotPos = [0.05 0.76 0.90 0.08];
        pplotPos = [0.05 0.64 0.90 0.08];
        vplotPos = [0.05 0.50 0.90 0.10];
        aplotPos = [0.05 0.36 0.90 0.10];
        fixplotPos = [0.05 0.04 0.43 0.28];
        rawplotPos = [0.52 0.04 0.43 0.28];
    else
        xplotPos = [0.05 0.88 0.90 0.08];
        yplotPos = [0.05 0.76 0.90 0.08];
        vplotPos = [0.05 0.60 0.90 0.12];
        aplotPos = [0.05 0.44 0.90 0.12];
        fixplotPos = [0.05 0.06 0.43 0.34];
        rawplotPos = [0.52 0.06 0.43 0.34];
    end
else
    if isfield(data,'pupil') && ~isempty(data.pupil.size)
        xplotPos = [0.05 0.88 0.90 0.08];
        yplotPos = [0.05 0.76 0.90 0.08];
        pplotPos = [0.05 0.60 0.90 0.12];
        vplotPos = [0.05 0.44 0.90 0.12];
        fixplotPos = [0.05 0.06 0.43 0.34];
        rawplotPos = [0.52 0.06 0.43 0.34];
    else
        xplotPos = [0.05 0.84 0.90 0.12];
        yplotPos = [0.05 0.68 0.90 0.12];
        vplotPos = [0.05 0.52 0.90 0.12];
        fixplotPos = [0.05 0.06 0.43 0.40];
        rawplotPos = [0.52 0.06 0.43 0.40];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% plot X trace with fixation markers
ax = axes('position',xplotPos);
plotWithMark(time,xdata,{'k-'},...                                      % data (y,x), style
             'time (ms) - fixations',xlbl,titel,...                     % x-axis label, y-axis label, axis title
             missFlag{:}, ...                                           % color part of trace that is missing
             blinkMarks{:}, ...                                         % blink markers (if any)
             fixmarks{:} ...                                            % fixation markers (if any)
            );
axis([mmt(1) mmt(2) rect(1) rect(3)]);
axis ij


%%% plot Y trace with fixation markers
ay = axes('position',yplotPos);
plotWithMark(time,ydata,{'k-'},...                                      % data (y,x), style
             'time (ms) - fixations',ylbl,'',...                        % x-axis label, y-axis label, axis title
             missFlag{:}, ...                                           % color part of trace that is missing
             blinkMarks{:}, ...                                         % blink markers (if any)
             fixmarks{:} ...                                            % fixation markers (if any)
            );
axis([mmt(1) mmt(2) rect(2) rect(4)]);


%%% plot pupil size trace with blink markers
if isfield(data,'pupil') && ~isempty(data.pupil.size)
    % pupil size:
    % determine axis size
    psr = max(data.pupil.size)-min(data.pupil.size);
    axisSize = [];
    if psr~=0
        axisSize = [mmt(1) mmt(2) min(data.pupil.size)-.03*psr max(data.pupil.size)+.03*psr];
    end
    ap = axes('position',pplotPos);
    hold on;
    plotWithMark(time,data.pupil.size,{'k-'},...                        % data (y,x), style
                 'time (ms) - blinks',plbl,'',...                       % x-axis label, y-axis label, axis title
                 missFlag{:}, ...                                       % color part of trace that is missing
                 blinkMarks{:} ...                                      % blink markers (if any)
                );
    if ~isempty(axisSize)
        axis(axisSize)
    end
    % change of pupil size:
    pvdat = abs(data.pupil.dsize);
    % determine axis size
    axisSize = [];
    if max(max(pvdat))~=0
        axisSize = [mmt(1) mmt(2) 0 max(pvdat)*1.03];
    end
    apv= axes('position',pplotPos);
    hold on;
    % line at 0
    plot([time(1) time(end)],[0 0],'b');
    hold on;
    plotWithMark(time,pvdat,{'k-'},...                                  % data (y,x), style
                 'time (ms) - blinks',pvlbl,'',...                      % x-axis label, y-axis label, axis title
                 missFlag{:}, ...                                       % color part of trace that is missing
                 blinkMarks{:} ...                                      % blink markers (if any)
                );
    if ~isempty(axisSize)
        axis(axisSize);
    end
    if isfield(data,'blink') && isfield(data.blink,'peakDSizeThreshold')
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
plotVel(time,vel{1},vlbl{1},'vel',datatype,...
    missFlag,sacon,sacoff,glismarks,blinkMarks,mmt,...
    qSaccadeTemplateRefinement,saccadePeakVelocityThreshold,saccadeOnsetVelocityThreshold,glissadeSearchSamples,saccadeOffsetVelocityThreshold);
avx = axes('position',vplotPos);
plotVel(time,vel{2},vlbl{2},'velX',datatype,...
    missFlag,sacon,sacoff,glismarks,blinkMarks,mmt,...
    qSaccadeTemplateRefinement,saccadePeakVelocityThreshold,saccadeOnsetVelocityThreshold,glissadeSearchSamples,saccadeOffsetVelocityThreshold);
avy = axes('position',vplotPos);
plotVel(time,vel{3},vlbl{3},'velY',datatype,...
    missFlag,sacon,sacoff,glismarks,blinkMarks,mmt,...
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

%%% either plot cross correlation output with saccade and glissade markers,
%%% or use the space to plot acceleration, or fuck it
aaxs = [];  % empty if no acceleration plots
ac   = [];  % empty if no cross correlation output plot
if qSaccadeTemplate
    % determine axis size
    axisSize = [mmt(1) mmt(2) 0 min(2.5,max(data.deg.velXCorr))];    % xcorr values above 2.5 seem to only occur due to noise
    ac = axes('position',cplotPos);
    hold on;
    % line at 0
    plot([time(1) time(end)],[0 0],'b');
    hold on;
    plotWithMark(time,data.deg.velXCorr,{'k-'},...                          % data (y,x), style
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
    axis(axisSize);
elseif qHaveAcceleration
    av2 = axes('position',aplotPos);
    plotAcc(time,acc{1},albl{1},'vel', missFlag,sacon,sacoff,glismarks,blinkMarks,mmt);
    avx = axes('position',aplotPos);
    plotAcc(time,acc{2},albl{2},'velX',missFlag,sacon,sacoff,glismarks,blinkMarks,mmt);
    avy = axes('position',aplotPos);
    plotAcc(time,acc{3},albl{3},'velY',missFlag,sacon,sacoff,glismarks,blinkMarks,mmt);
    aaxs = [av2 avx avy];
    % show desired vel at start
    toHide = [1:3]; toHide(toHide==vidx) = [];
    for p=toHide
        set([aaxs(p); allchild(aaxs(p))],'visible','off');
    end
    axes(aaxs(vidx));   % set visible axis to current and topmost axis
end

% link x-axis (time) of the three/four timeseries for easy viewing
linkaxes([ax ay ap apv vaxs aaxs ac],'x');


%%% plot scanpath of raw data and of fixations
if qHaveFixations
    asf = axes('position',fixplotPos);
    hold on
    if nargin>=8 && strcmp(datatype,'pix') && ~isempty(pic)
        imagesc([0 size(pic.imdata,2)]+pic.offset(2),[0 size(pic.imdata,1)]+pic.offset(1),pic.imdata);
    end
    plotWithMark(xfixpos,yfixpos,{'k-'},...                                             % data (y,x), style
                 xlbl,ylbl,'',...                                                       % x-axis label, y-axis label, axis title
                 [1:length(xfixpos)],{'go','MarkerFaceColor','g','MarkerSize',4},...    % mark each fixation (that is marker on each datapoint we feed it
                 1,                  {'co','MarkerFaceColor','c','MarkerSize',4},...    % make first fixation marker blue
                 length(xfixpos),    {'mo','MarkerFaceColor','m','MarkerSize',4} ...    % make last  fixation marker red
                );
    axis(rect([1 3 2 4]));
else
    asf = [];
end

asr = axes('position',rawplotPos);
hold on
if nargin>=8 && strcmp(datatype,'pix') && ~isempty(pic)
    imagesc([0 size(pic.imdata,2)]+pic.offset(2),[0 size(pic.imdata,1)]+pic.offset(1),pic.imdata);
end
plotWithMark(xdata,ydata,{'k-'},...                                                 %  data (y,x), style
             xlbl,ylbl,'',...                                                       % x-axis label, y-axis label, axis title
             1,                  {'co','MarkerFaceColor','c','MarkerSize',4},...    % use blue marker for first datapoint
             length(xdata),      {'mo','MarkerFaceColor','m','MarkerSize',4} ...    % use red  marker for last  datapoint
            );
axis(rect([1 3 2 4]));

% link view of the two scanpath plots for easy viewing
linkaxes([asr asf],'xy');

% make sure we don't lose the standard toolbar
set(gcf,'Toolbar','figure');
set(gcf,'DockControls','off');
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
        qShow   = strcmp(but,strs);
        toHidev = find(~qShow);
        toShow  = find( qShow);
        % do show and hide
        for q=toHidev
            set([vaxs(q); allchild(vaxs(q))],'visible','off');
            if ~isempty(aaxs)
                set([aaxs(q); allchild(aaxs(q))],'visible','off');
            end
        end
        set([vaxs(toShow); allchild(vaxs(toShow))],'visible','on');
        axes(vaxs(toShow));     % set visible axis to current and topmost axis
        if ~isempty(aaxs)
            set([aaxs(toShow); allchild(aaxs(toShow))],'visible','on');
            axes(aaxs(toShow));     % set visible axis to current and topmost axis
        end
        % adjust labels on buttons
        strs2v  = strs(toHidev);
        set(vt1,'String',strs2v{1});
        set(vt2,'String',strs2v{2});
    end

end

function plotVel(time,vel,vlbl,veltype,datatype,...
    missFlag,sacon,sacoff,glisMarks,blinkMarks,mmt,...
    qSaccadeTemplateRefinement,saccadePeakVelocityThreshold,saccadeOnsetVelocityThreshold,glissadeSearchSamples,saccadeOffsetVelocityThreshold)
% determine axis size
axisSize = [];
if any(~isnan(vel))
    if min(0,min(vel))==0
        axisSize = [mmt(1) mmt(2) 0 max(vel)*1.03];
    else
        psr = max(vel)-min(vel);
        axisSize = [mmt(1) mmt(2) min(vel)-.03*psr max(vel)+.03*psr];
    end
end
hold on;
% line at 0
plot([time(1) time(end)],[0 0],'b');
plotWithMark(time,vel,{'k-'},...                                        % data (y,x), style
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
if ~isempty(axisSize)
    axis(axisSize)
end
if ~strcmp(veltype,'vel')
    axis ij
end
end

function plotAcc(time,acc,albl,veltype,...
    missFlag,sacon,sacoff,glisMarks,blinkMarks,mmt)
% determine axis size
axisSize = [];
if any(~isnan(acc))
    if min(0,min(acc))==0
        axisSize = [mmt(1) mmt(2) 0 max(acc)*1.03];
    else
        psr = max(acc)-min(acc);
        axisSize = [mmt(1) mmt(2) min(acc)-.03*psr max(acc)+.03*psr];
    end
end
hold on;
% line at 0
plot([time(1) time(end)],[0 0],'b');
plotWithMark(time,acc,{'k-'},...                                        % data (y,x), style
             'time (ms) - saccades/glissades',albl,'',...               % x-axis label, y-axis label, axis title
             missFlag{:}, ...                                           % color part of trace that is missing
             sacon, {'bo','MarkerFaceColor','blue','MarkerSize',4},...  % saccade on  markers
             sacoff,{'ro','MarkerFaceColor','red' ,'MarkerSize',4},...  % saccade off markers
             glisMarks{:}, ...                                          % glissade markers (if any)
             blinkMarks{:} ...                                          % blink markers (if any)
    );
if ~isempty(axisSize)
    axis(axisSize)
end
if ~strcmp(veltype,'vel')
    axis ij
end
end
