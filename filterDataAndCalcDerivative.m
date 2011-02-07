function data = filterDataAndCalcDerivative(data,ETparams)

% prepare parameters
%--------------------------------------------------------------------------
% span of filter, use minimum length of saccade
span = ceil(ETparams.saccade.minDur/1000*ETparams.samplingFreq);
% number of tabs of filter
F = 2*ceil(span)-1;

% compute eye positions, velocities, and accelerations in pixels, convert
% to degree
%--------------------------------------------------------------------------

% Calculate the filtered position, velocity and acceleration
[tempP,tempV,tempA] = sgFilt([data.deg.Xori data.deg.Yori],[0 1 2],F);

% store and calculate derivative magnitudes
data.deg.X      = tempP(:,1);
data.deg.Y      = tempP(:,2);

data.deg.vel    = hypot(tempV(:,1), tempV(:,2)) * ETparams.samplingFreq;
data.deg.acc    = hypot(tempA(:,1), tempA(:,2)) * ETparams.samplingFreq^2;




%%% helpers
function varargout = sgFilt(x,difforder,ntaps)
% wrapper for convenient syntax, fitting always using 2nd order polynomial
for p=1:length(difforder)
    varargout{p} = savitzkyGolayFilt(x,2,difforder(p),ntaps);
end