function data = prepareData(x,y,pupilsize,ETparams)
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
% positive X axis points rightward and the positive Y axis points downward.
% Furthermore, positive rotations around the Y axis are rightward, but
% positive rotations around the X axis are downward too (frankly, this is
% not a consistent right-handed or left-handed sytem, this rotation should
% then be upward). This doesn't matter until we want to correctly calculate
% the axis of eye rotation. We'll then have to revise all code to conform
% to Haslwanter's convention, for ease of use.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Positive axes with Z pointing from head into the screen
%
%                 Z
%                #
%              ###
%            #####
%            #
%           #
%          #
%        #
%       #                                     #
%      #                                      ###
%     ############################################ X
%     #                                       ###
%     #                                       #
%     #
%     #
%     #
%     #
%     #
%     #
%     #
%     #
%     #
%     #
%   #####
%    ###
%     #
%     Y

% ensure column vectors
data.pix.X      = x(:);
data.pix.Y      = y(:);
data.pupil.size = pupilsize(:);

% flip X if specified
if ETparams.data.qFlipX
    data.pix.X  = -data.pix.X + ETparams.screen.resolution(1);
end
% flip Y if specified
if ETparams.data.qFlipY
    data.pix.Y  = -data.pix.Y + ETparams.screen.resolution(2);
end

% move origin to point on screen straight ahead of subject
data.pix.X      = data.pix.X - ETparams.screen.subjectStraightAhead(1);
data.pix.Y      = data.pix.Y - ETparams.screen.subjectStraightAhead(2);

% finally, convert gaze position in pixels to Fick angles in degree
[data.deg.Azi, data.deg.Ele] = pix2fick(data.pix.X,data.pix.Y,ETparams);