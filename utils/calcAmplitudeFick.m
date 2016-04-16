function [amplitude,direction] = calcAmplitudeFick(h1,v1,h2,v2)

% This function calculates the angular distance between two points that are
% specified in Fick angles. It provides an exact way of computing the
% angular distance between two gaze points, using the dot product of the
% two gaze vectors to calculate it.
% input:
%   h1: Fick azimuth   angle of direction 1, in degree
%   v1: Fick elevation angle of direction 1, in degree
%   h2: Fick azimuth   angle of direction 2, in degree
%   v2: Fick elevation angle of direction 2, in degree
% output:
%   amplitude: angle (in degree) between the two gaze directions specified
%              by (h1,v1) and (h2,v2)
%   direction: orientation (in degree) of eye movement on the tangent plane
%              (on the screen)


d2r = pi/180.0;
% get unit vector in direction of gaze points (note labeling of
% outputs, in my mind x is leftward and y is upward)
% note: MATLAB's sph2cart and cart2sph use Fick angle order of rotations
[z1,x1,y1] = sph2cart(h1(:)*d2r, v1(:)*d2r, 1);
[z2,x2,y2] = sph2cart(h2(:)*d2r, v2(:)*d2r, 1);

% compute angle between them using dot product (simple as we already
% know vectors have unit length)
amplitude = acosd(dot([x1,y1,z1],[x2,y2,z2],2));

if nargout==2
    % compute direction on tangent screen
    % (divide by z to project to screen at 1m)
    % negate y coordinate as on a screen, increasing y is downward
    direction = atan2(-(y2./z2-y1./z1),x2./z2-x1./z1) * 180/pi;
end