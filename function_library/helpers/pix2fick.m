function [fx,fy,dist] = pix2fick(px,py,varargin)

r2d = 180.0/pi;

% convert gaze position in pixels to Fick angles in degree
% first convert pixels to cm away from origin
if nargin==3
    pixPerMeter     = varargin{1}.screen.resolution ./ varargin{1}.screen.size;
    viewDist        = varargin{1}.screen.viewingDist;
else
    pixPerMeter     = varargin{1} ./ varargin{2}*100;   % resolution/size
    viewDist        = varargin{3}/100;                  % viewing distance
end

% Then convert to Fick angles (MATLAB's cart2sph models a Fick gimbal)
% although their reference Z axis is our Y axis, their X axis is our Z axis
% and their Y axis is our X axis:
% cart2sph: X Y Z
% our Fick: Z X Y
[fx,fy,dist]    = cart2sph(viewDist, px./pixPerMeter(1), py./pixPerMeter(2));

% convert to degrees (Fick angles)
fx              = fx*r2d;
fy              = fy*r2d;
