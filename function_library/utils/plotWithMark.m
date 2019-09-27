function hndls = plotWithMark(xdata,ydata,addInps,usrData,xlbl,ylbl,titel,varargin)

% makes and x-y plot with specified axis labels and title in the current
% axes. addInps specifies additional input to the plot command that draws the x-y data.
% Optionally plots markers, specified by the varargin part
%
% These markers are identified by pairs of input. You can specify an
% unlimited number of those pairs. The first is indices into the xdata and
% ydata, the second is a cell containing 0-x additional inputs to the
% plot() command, e.g. setting the look of the markers
%
% plotWithMark(...,markerIdxs,{});
% plotWithMark(...,markerIdxs,{'ro'});
% plotWithMark(...,markerIdxs,{'bo','MarkerFaceColor','blue','MarkerSize',4};
%
% are thus all valid calls.

narginchk(7, inf);
% length of varargin must be even as marker inputs come in pairs
assert(mod(length(varargin),2)==0,'Number of inputs related to markers must be even as they come in pairs, got %d inputs related to markers',nargin-5);
if isempty(addInps)
    addInps = {'k-'};
end
usrData.x = xdata;
usrData.y = ydata;

% plot trace
hndls = plot(xdata,ydata,addInps{:},'LineWidth',1,'UserData',usrData);
if ~isempty(xlbl)
    xlabel(xlbl);
end
if ~isempty(ylbl)
    ylabel(ylbl);
end
if ~isempty(titel)
    title(titel);
end

if isempty(xdata)
    % nothing will be plotted anyway and indexing below will crash
    return;
end

if isstruct(usrData)
    fields = fieldnames(usrData);
    fields(strcmp(fields,'tag')) = [];
end

hold on;
for p=1:2:length(varargin)
    usrDatap = usrData;
    qUserDataNoGrow = strcmp(varargin{p+1},'UserDataNoGrow');
    if any(qUserDataNoGrow)
        idxs = find(qUserDataNoGrow);
        qUserDataNoGrow = ~~varargin{p+1}{idxs(end)+1};
        if qUserDataNoGrow
            usrDatap.dontGrow = true;
        end
        varargin{p+1}([idxs; idxs+1]) = [];
    end
    if isstruct(usrData)
        for f=1:length(fields)
            usrDatap.(fields{f}) = usrDatap.(fields{f})(varargin{p});
        end
    end
    h = plot(xdata(varargin{p}),ydata(varargin{p}),varargin{p+1}{:},'UserData',usrDatap);
    if ~isempty(h)
        hndls = [hndls h]; %#ok<AGROW>
    end
end
hold off;
