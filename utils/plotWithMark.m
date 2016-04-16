function plotWithMark(xdata,ydata,addInps,xlbl,ylbl,titel,varargin)

% makes and x-y plot with specified axis lalbels and title in the current
% axes. Optionally plots markers, specified by the varargin plot
%
% Markers are identified by pairs of input. You can specify an unlimited
% number of those pairs. The first is indices into the xdata and ydata, the
% second is a cell containing 0-x additional inputs to the plot() command
%
% plotWithMark(...,marker,{});
% plotWithMark(...,marker,{'ro'});
% plotWithMark(...,marker,{'bo','MarkerFaceColor','blue','MarkerSize',4};
%
% are thus all valid calls.

narginchk(6, inf);
% length of varargin must be even as marker inputs come in pairs
assert(mod(length(varargin),2)==0,'Number of inputs related to markers must be even as they come in pairs, got %d inputs related to markers',nargin-5);
if isempty(addInps)
    addInps = {'k-'};
end

% plot trace
plot(xdata,ydata,addInps{:},'LineWidth',1);
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

hold on;
for p=1:2:length(varargin)
    plot(xdata(varargin{p}),ydata(varargin{p}),varargin{p+1}{:});
end
hold off;