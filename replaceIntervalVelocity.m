function [vel,velX,velY] = replaceIntervalVelocity(vel,velX,velY,on,off)
% Linearly interpolate velocity between on and off
% on and off are sample numbers (data indices) of beginning and end of
% interval to be replaced by linear interpolation between begin and end

npoint = off-on+1;

% components: replace with linearly interpolated velocity
velX(on:off) = linspace(velX(on),velX(off),npoint);
velY(on:off) = linspace(velY(on),velY(off),npoint);

% calculate interpolate 2D velocity. In effect this is now interpolated
% with a bicubic spline. Thats fine, good even as no edges are introduced
% into the data, as we care most about the component velocities in the
% situations I can think of.
vel(on:off)  = hypot(velX(on:off),velY(on:off));