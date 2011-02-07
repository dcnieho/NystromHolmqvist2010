function data = prepareData(x,y,ETparams)
% prepares data for the rest of the algorithm. Moves origin, converts to
% degrees and things like that.

% ensure column vectors
data.pix.X = x(:);
data.pix.Y = y(:);

% TODO: move origin to center of screen, unless option says something more
% complicated

% TODO flip Y if specified

% TODO convert to degrees
% if assumeTangentLinearity
% Calculate how many degrees one pixel spans.
[angleInPixelsH, angleInPixelsV] = degrees2pixels(1, ETparams.screen.viewingDist, ETparams.screen.resolution, ETparams.screen.size);
data.deg.Xori   = data.pix.X / angleInPixelsH;
data.deg.Yori   = data.pix.Y / angleInPixelsV;
% else complicated model to compute eccentricity from straight ahead
% correctly using cosine rule