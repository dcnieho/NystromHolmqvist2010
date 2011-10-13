function [px,py] = fick2pix(fx,fy,ETparams)

d2r = pi/180.0;
% first convert pixels to cm away from origin
pixPerMeter = ETparams.screen.resolution ./ ETparams.screen.size;

% create Cartesian unit vector
[z,px,py]   = sph2cart(fx*d2r, fy*d2r, 1);

% project to screen at screen's distance
px          = px./z.*ETparams.screen.viewingDist;
py          = py./z.*ETparams.screen.viewingDist;

% convert cm away from origin to pixels
px          = px*pixPerMeter(1);
py          = py*pixPerMeter(2);