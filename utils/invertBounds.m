function [ion,ioff] = invertBounds(on,off,max)
% return bounds for the stretches not enclosed by the input bounds

ion = off(1:end-1)+1;
ioff= on(2:end)-1;

if on(1)~=1
    ion = [1 ion];
    ioff = [on(1)-1 ioff];
end
if off(end)~=max
    ion = [ion off(end)+1];
    ioff= [ioff max];
end