function [fx,fy,dist] = pix2fick(px,py,ETparams)

r2d = 180.0/pi;

% convert gaze position in pixels to Fick angles in degree
% first convert pixels to cm away from origin
pixPerMeter     = ETparams.screen.resolution ./ ETparams.screen.size;

% Then convert to Fick angles (MATLAB's cart2sph models a Fick gimbal)
% although their reference Z axis is our Y axis, their X axis is our Z axis
% and their Y axis is our X axis:
% cart2sph: X Y Z
% our Fick: Z X Y
[fx,fy,dist]    = cart2sph(ETparams.screen.viewingDist, px./pixPerMeter(1), py./pixPerMeter(2));

% convert to degrees (Fick angles)
fx              = fx*r2d;
fy              = fy*r2d;