function data = prepareData(x,y,ETparams)
% prepares data for the rest of the algorithm. Moves origin, converts to
% degrees and things like that.

% ensure column vectors
data.pix.Xori = x(:);
data.pix.Yori = y(:);

% flip X if specified
if ETparams.data.qFlipX
    data.pix.Xori = -data.pix.Xori + ETparams.screen.resolution(1);
end
% flip Y if specified
if ETparams.data.qFlipY
    data.pix.Yori = -data.pix.Yori + ETparams.screen.resolution(2);
end

% move origin to point on screen straight ahead of subject
data.pix.Xori = data.pix.Xori - ETparams.screen.subjectStraightAhead(1);
data.pix.Yori = data.pix.Yori - ETparams.screen.subjectStraightAhead(2);