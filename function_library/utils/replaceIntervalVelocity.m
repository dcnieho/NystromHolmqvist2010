function [vel,velX,velY] = replaceIntervalVelocity(vel,velX,velY,Y,qDeg,on,off,val)
% Linearly interpolate velocity between on and off
% on and off are sample numbers (data indices) of beginning and end of
% interval to be replaced by linear interpolation between begin and end, or
% by constant value if input val is set
% if interval to replace ends at end, we can't do anything, unless fixed
% value requested
if off==length(vel) && nargin<8
    velX(on:off) = nan;
    velY(on:off) = nan;
     vel(on:off) = nan;
    return;
end


npoint = off-on+1;

if nargin>7
    % components: replace with constant velocity
    velX(on:off) = val;
    velY(on:off) = val;
    if qDeg
        vel(on:off)  = sqrt(val.^2.*cosd(Y(on:off)) + val.^2);
    else
        vel(on:off)  = hypot(val,val);
    end
    
else
    % components: replace with linearly interpolated velocity
    velX(on:off) = linspace(velX(on),velX(off),npoint);
    velY(on:off) = linspace(velY(on),velY(off),npoint);
    
    % calculate interpolate 2D velocity. In effect this is now interpolated
    % with a bicubic spline. Thats fine, good even as no edges are introduced
    % into the data, as we care most about the component velocities in the
    % situations I can think of.
    if qDeg
        vel(on:off)  = sqrt(velX(on:off).^2.*cosd(Y(on:off)) + velY(on:off).^2);
    else
        vel(on:off)  = hypot(velX(on:off),velY(on:off));
    end
end
