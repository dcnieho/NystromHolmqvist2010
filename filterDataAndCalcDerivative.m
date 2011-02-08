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
%
% Eye position data is stored in Fick (1854) angles. Given a
% head-referenced frame in which the Z axis points upwards (dorso-ventral),
% the X axis forward (naso-occipital) and the Y axis leftward
% (inter-aural), Fick angles are given by a rotation around Z axis first
% (yaw/azimuth rotation), followed by rotation around the Y axis
% (pitch/elevation rotation) and lastly a rotation around the X axis
% (roll/torsion rotation).
% Eye rotation velocity and acceleration are currently computed assuming a
% linear relationship between position on the tangent plane and angle,
% simply as I don't know how to compute them correctly in Fick angles. This
% is accurate if gaze is pointed towards straight ahead, but it means
% velocity and acceleration for eccentric eye positions are overestimated.
% Fick A (1854). Die bewegungen des menschlichen augapfels. Zeitschrift für
% rationelle Medizin 4: 109-128. 

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

% Calculate the filtered position, velocity and acceleration
[tempP,tempV,tempA] = sgFilt([data.pix.Xori data.pix.Yori],[0 1 2],ntaps);

% store filtered
data.pix.X      = tempP(:,1);
data.pix.Y      = tempP(:,2);

% and convert to Fick angles in degree
% first convert pixels to cm away from origin
pixPerMeter = ETparams.screen.resolution ./ ETparams.screen.size;
% Then convert to Fick angles (MATLAB's cart2sph happens to use that order)
[dx,dy]     = cart2sph(ETparams.screen.viewingDist,data.pix.X ./ pixPerMeter(1),data.pix.Y ./ pixPerMeter(2));
% convert to degrees (Fick angles)
data.deg.X  = dx./pi*180;
data.deg.Y  = dy./pi*180;

% calculate derivative magnitudes
if ETparams.data.qPreciseCalcDeriv
    % TODO: implement
    % how to compute eye velocity and acceleration in Fick angles using the
    % output from the derivative filter?
else
    % calculate pixels per degree
    meterPerDegree  = ETparams.screen.viewingDist * tand(1);    % 1° movement away from straight ahead
    pixPerDegree    = pixPerMeter * meterPerDegree;             % dimensional analysis of the variable names show this is correct :P
    % get magnitude of velocity vector in pixels and convert to degree
    % assuming linearity (good enough for most purposes / small screens)
    % dont forget to scale by sampling rate.
    data.deg.vel    = hypot(tempV(:,1)/pixPerDegree(1), tempV(:,2)/pixPerDegree(2)) * ETparams.samplingFreq;
    data.deg.acc    = hypot(tempA(:,1)/pixPerDegree(1), tempA(:,2)/pixPerDegree(2)) * ETparams.samplingFreq^2;
end



%%% helpers
function varargout = sgFilt(x,difforder,ntaps)
% wrapper for convenient syntax, fitting always using 2nd order polynomial
for p=1:length(difforder)
    varargout{p} = savitzkyGolayFilt(x,2,difforder(p),ntaps);
end