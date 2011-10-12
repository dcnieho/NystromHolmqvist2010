function data = prepareData(x,y,ETparams)
% prepares data for the rest of the algorithm. Moves origin, flip data,
% converts to Fick angles in degree, things like that.
%
% Eye position data is stored in Fick (1854) angles. Given a
% head-referenced frame in which the Y axis points upwards (dorso-ventral),
% the Z axis forward (naso-occipital) and the X axis leftward
% (inter-aural), Fick angles are given by a rotation around Y axis first
% (yaw/azimuth rotation), followed by rotation around the X axis
% (pitch/elevation rotation) and lastly a rotation around the Z axis
% (roll/torsion rotation).
% Fick A (1854). Die bewegungen des menschlichen augapfels. Zeitschrift für
% rationelle Medizin 4: 109-128.
% See also:
% Haslwanter T (1995) Mathematics of 3-dimensional eye rotations. Vision
% Res 35, 1727-1739, e.g. Figure 4
% 
% Furthermore, it should be noted that due to our lab conventions, the
% positive X axis points righward and the positive Y axis points downward.
% Furthermore, positive rotations around the Y axis are rightward, but
% positive rotations around the X axis are downard too (frankly, this is
% not a consistent right-handed sytem, this rotation should then be
% upward). This doesn't matter until we want to correctly calculate the teh
% axis of eye rotation. We'll then have to revise all code to conform to
% Haslwanter's convention, for ease of use.

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

% finally, convert gaze position in pixels to Fick angles in degree
% first convert pixels to cm away from origin
pixPerMeter = ETparams.screen.resolution ./ ETparams.screen.size;
% Then convert to Fick angles (MATLAB's cart2sph models a Fick gimbal)
% although their reference Z axis is our Y axis, their X axis is our Z axis
% and their Y axis is our X axis:
% cart2sph: X Y Z
% our Fick: Z X Y
[dx,dy]     = cart2sph(ETparams.screen.viewingDist, data.pix.X./pixPerMeter(1), data.pix.Y./pixPerMeter(2));
% convert to degrees (Fick angles)
data.deg.X  = dx./pi*180;
data.deg.Y  = dy./pi*180;