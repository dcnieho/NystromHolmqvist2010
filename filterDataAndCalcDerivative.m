function data = filterDataAndCalcDerivative(data,ETparams)
% We'll use the Savitzky-Golay filter for smoothing and differentiation,
% which basically performs a polynomial least square fit in a moving
% window. The eye position data is smoothed by extracting the zeroth order
% coefficient of this polynomial. Eye velocity and acceleration and
% obtained by returning the higher order coefficients of the fitted
% polynomial. Note that for kth order derivatives (k>1) we'd need to
% multiply the ouput by k! to obtain the derivative coefficient (see Taylor
% series expansion). This already taken into account internally in the
% savitzkyGolay function (we use its B output, not its G output for which
% such multiplication would be needed), so you don't need to worry about
% it. We'll need to divide the outputted value by T^K (where T is the
% sampling interval as everything is a function of time) though to get the
% correct value of the derivative.
% See also sections 5.7 (on calculating numerical derivatives) and 14.8 (on
% Savitzky-Golay filters) in Numerical recipes in C++
%
% Note that Numerical Recipes recommends using a polynomial order of at
% least 4 when you're interested in derivatives of the data. Savitzky-Golay
% filters however introduce false negative peaks at the edges of large
% positive peaks. These might then be detected as false glissades in our
% algorithm as they become positive peaks after the saccade when
% calculating the 2D velocity magnitude. We therefore use a second order
% polynomial only, which has LESS of this problem, but is still affected by
% it. I am not sure about the extend to which this increases the number of
% glissades found, i.e., how many of the detected glissades are just
% artifacts of the filter. We might look at weighting the polynomial
% regression filter coefficients, as supported by the savitzkyGolay()
% function.
% This brings with it another problem. Savitzki-Golay filters of polynomial
% order N preserve up to the Nth order moment of their input data. 0th
% moment: area under a peak, 1st order: position of the peak, 2nd order
% moment: width of the peak. I can't find more, but I'd bet that the 3th
% order is the width of the velocity peak. This means that although the
% location of the velocity peak is preserved, its width is not. Indeed,
% comparing a first order polynomial filter on position to the performance
% of a second order filter, we see a similar change in the shape of peaks
% of the position signal as when comparing the shape of the velocity
% signal's peaks using a second and a third order filter. This means that,
% in addition to possibly introducing extra glissades into the data, we
% also catch slightly longer saccades of slightly lower velocity using the
% Savitzky-Golay filter.
% Other filtering methods might bring resolve here, such as possibly an
% optimal (Wiener) filter, or some other more advanced and recently
% developed methods. I am not knowledgable on this subject matter though,
% so I'll need to find a DSP guy to resolve this with once and for all
% (hopefully).

% prepare parameters
%--------------------------------------------------------------------------
% span of filter, use minimum length of saccade. Its very important to not
% make the filter window wider than the narrowest feature we are interested
% in, or we'll smooth out those features.
window  = ceil(ETparams.saccade.minDur/1000*ETparams.samplingFreq);
% number of filter taps
ntaps   = 2*ceil(window)-1;

% compute eye positions, velocities, and accelerations in pixels, convert
% to degree
%--------------------------------------------------------------------------

if ETparams.data.qPreciseCalcDeriv
    % TODO implement: unfinished below!
    % Calculate the filtered position
    tempP = sgFilt([data.deg.Xori data.deg.Yori],0,ntaps);
    
    % store filtered
    data.deg.X      = tempP(:,1);
    data.deg.Y      = tempP(:,2);
    
    % Now calculate eye velocity and acceleration precisely
    % Create quaternion rotation vector denoting eye position for each
    % sample, then take derivative of the series of these. (Note: this does
    % not use the analytical derivatives of the polynomial filtering to
    % compute eye velocity/acceleration (TODO: derive and see if we can do
    % it!), but uses numerical derivatives of our filtered eye position.
    % See e.g.:
    % 
else
    % Calculate the filtered position, velocity and acceleration
    [tempP,tempV,tempA] = sgFilt([data.deg.Xori data.deg.Yori],[0 1 2],ntaps);
    
    % store filtered
    data.deg.X      = tempP(:,1);
    data.deg.Y      = tempP(:,2);
    
    % calculate derivative magnitudes
    % note NOTE NOTE: This works fine for our purpose of detecting velocity
    % and acceleration peaks, which in essense is all the algorithm really
    % does. Both fixations and smooth pursuit are periods of low velocity
    % interspersed with high velocity peaks of the saccades. So as long as
    % we get good peaks we can run our algorithm and we don't have to care
    % about their exact height. However, when you are interested in the
    % exact eye velocity/acceleration (both during saccades and during
    % pursuit), you'll have to use the derivatives of the eye rotation
    % vectors (the exact method above). This trick using Pythagoras'
    % theorem is pretty crude as actually it only applies in Cartesian
    % space, not in the spherical system we use.
    % get magnitude of velocity and acceleration vectors, dont forget to
    % scale by sampling rate.
    data.deg.vel    = hypot(tempV(:,1), tempV(:,2)) * ETparams.samplingFreq;
    data.deg.acc    = hypot(tempA(:,1), tempA(:,2)) * ETparams.samplingFreq^2;
end



%%% helpers
function varargout = sgFilt(x,difforder,ntaps)
% wrapper for convenient syntax, fitting always using 2nd order polynomial
for p=1:length(difforder)
    varargout{p} = savitzkyGolayFilt(x,2,difforder(p),ntaps);
end