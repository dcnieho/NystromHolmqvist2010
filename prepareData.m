function data = prepareData(x,y,ETparams)
% prepares data for the rest of the algorithm. Moves origin, converts to
% degrees and things like that.

% ensure column vectors
data.pix.X = x(:);
data.pix.Y = y(:);

% flip X if specified
if ETparams.data.qFlipX
    data.pix.X = -data.pix.X + ETparams.screen.resolution(1);
end
% flip Y if specified
if ETparams.data.qFlipY
    data.pix.Y = -data.pix.Y + ETparams.screen.resolution(2);
end

% move origin to point on screen straight ahead of subject
data.pix.X = data.pix.X - ETparams.screen.subjectStraightAhead(1);
data.pix.Y = data.pix.Y - ETparams.screen.subjectStraightAhead(2);

% convert eye position to eccentricity from straight ahead in degrees
% convert pixels to cm away from origin
pixPerMeter = ETparams.screen.resolution ./ ETparams.screen.size;
px          = data.pix.X ./ pixPerMeter(1);
py          = data.pix.Y ./ pixPerMeter(2);
% convert to spherical system in radian
[dx,dy] = cart2sph(ETparams.screen.viewingDist,px,py);

% convert to degrees
data.deg.Xori   = dx./pi*180;
data.deg.Yori   = dy./pi*180;