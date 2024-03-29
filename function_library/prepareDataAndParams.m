function [data,ETparams] = prepareDataAndParams(x,y,pupilsize,ETparams,time)
% prepares data for the rest of the algorithm. Moves origin, flip data,
% converts to Fick angles in degree, things like that.
% also deals with tracker's particular method of indicating missing data
% (e.g. outputting (0,0) and replacing that with nan)
%
% Eye position data is stored in Fick (1854) angles. Given a
% head-referenced frame in which the Y axis points upwards (dorso-ventral),
% the Z axis forward (naso-occipital) and the X axis leftward
% (inter-aural), Fick angles are given by a rotation around Y axis first
% (yaw/azimuth rotation), followed by rotation around the X axis
% (pitch/elevation rotation) and lastly a rotation around the Z axis
% (roll/torsion rotation).
% Fick A (1854). Die bewegungen des menschlichen augapfels. Zeitschrift f�r
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

if nargin>4
    data.time       = time(:);
end

% first flag missing data for this tracker
if ~isempty(ETparams.data.missCoords)
    qMissing = data.pix.X==ETparams.data.missCoords(1) & data.pix.Y==ETparams.data.missCoords(2);
    
    if any(qMissing)
        % now remove data from all eye position fields that we have
        data.pix    = replaceElementsInStruct(data.pix,qMissing,nan);
    end
end

% translate tracker output into screen coordinates
data.pix.X      = data.pix.X + ETparams.screen.resolution(1)/2 - ETparams.screen.dataCenter(1);
data.pix.Y      = data.pix.Y + ETparams.screen.resolution(2)/2 - ETparams.screen.dataCenter(2);

% flip X if specified
if ETparams.data.flipX
    data.pix.X  = -data.pix.X + ETparams.screen.resolution(1);
end
% flip Y if specified
if ETparams.data.flipY
    data.pix.Y  = -data.pix.Y + ETparams.screen.resolution(2);
end

% move origin to point on screen straight ahead of subject
data.pix.X      = data.pix.X - ETparams.screen.subjectStraightAhead(1);
data.pix.Y      = data.pix.Y - ETparams.screen.subjectStraightAhead(2);

% finally, convert gaze position in pixels to Fick angles in degree
[data.deg.Azi, data.deg.Ele, data.deg.dist] = pix2fick(data.pix.X,data.pix.Y,ETparams);

%%%% params
if isempty(data.pupil.size)
    % make sure classification of blinks by pupil size is not requested if
    % we don't have pupil size data
    ETparams.blink.classifyMode = ETparams.blink.classifyMode - bitand(ETparams.blink.classifyMode,uint8(1));
end
