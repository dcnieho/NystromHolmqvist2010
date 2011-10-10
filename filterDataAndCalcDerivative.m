function data = filterDataAndCalcDerivative(data,ETparams)
% By default, we'll use the Savitzky-Golay filter for smoothing and
% differentiation (simply numerical difference calculation is also
% possible), which basically performs a polynomial least square fit in a
% moving window. The eye position data is smoothed by extracting the zeroth
% order coefficient of this polynomial. Eye velocity and acceleration and
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
% In general, we'll need to asses how problematic this really is. When it
% comes to the filter, we need to simulate its performance (how much
% widening, how much undershoot is there really?). Generate a function much
% like the eye movement, but without any noise (Gaussian velocity profile,
% very low std for sharp peak? Or more directly cut a part from some
% eyemovement data and smooth it by hand) and see what happens to it when
% filtering by polynomial fit... Also, remember the Luo 2005 filter on
% properties of Savitzky-Golay differentiation filters. By the way, i
% indeed think my logic above was correct: the zeroth order moment of the
% first derivative of a time series should correspond/be relatable to the
% first order moment of the time series itself. It is thus a problem that
% the width of the velocity peak (3rd moment) is not present by a second
% order polynomial fit.

% calculate component velocities and accelerations
if ETparams.data.qNumericallyDifferentiate
    tempV   = diff([data.deg.X data.deg.Y],1,1);
    tempA   = diff([data.deg.X data.deg.Y],2,1);
    
    % make same length as position trace
    tempV   = [tempV; tempV( end,:)];
    tempA   = [tempA; tempA([end end],:)];
    
    % also calculate derivatives for eye position in pixels
    if ETparams.data.qAlsoStoreandSmoothPixels
        tempVpix    = diff([data.pix.X data.pix.Y],1,1);
        tempApix    = diff([data.pix.X data.pix.Y],2,1);
        
        % make same length as position trace
        tempVpix    = [tempVpix; tempVpix( end,:)];
        tempApix    = [tempApix; tempApix([end end],:)];
    end
else
    % prepare parameters
    %--------------------------------------------------------------------------
    % span of filter, use minimum length of saccade. Its very important to not
    % make the filter window wider than the narrowest feature we are interested
    % in, or we'll smooth out those features.
    window  = ceil(ETparams.saccade.minDur/1000*ETparams.samplingFreq);
    % number of filter taps
    ntaps   = 2*ceil(window)-1;
    
    % calculate derivatives
    [tempV,tempA] = sgFilt([data.deg.X data.deg.Y],[1 2],ntaps);
    
    % also calculate derivatives for eye position in pixels
    if ETparams.data.qAlsoStoreandSmoothPixels
        [tempVpix,tempApix] = sgFilt([data.pix.X data.pix.Y],[1 2],ntaps);
    end
end

% compute eye velocities, and accelerations
%--------------------------------------------------------------------------
if ETparams.data.qPreciseCalcDeriv
    % TODO: implement
    error('TODO')
    
    % Now calculate eye velocity and acceleration precisely
    % See Equation 30 in Haslwanter T (1995) Mathematics of 3-dimensional
    % eye rotations, Vision Res 35, 1727-1739. (look up that original
    % Goldstein reference as well). Might also want to look up Fetter M,
    % Haslwanter T, Misslisch M, Tweed D (1997) Three-dimensional kinematic
    % principles of eye-, head-, and limb movements, Harwood Academic
    % Publishers: Amsterdam, sounds relevant.
    % This allows us to use the analytical differentiation of the
    % polynomial fitted to eye azimuth and elevation (so-called coordinate
    % velocities). This however still doesn't give me a scalar velocity I
    % think.. need to look at the details of this and familiarize myself
    % with the mathematical tools. Once we got velocity, I think we can
    % safely take the acceleration numerically without adding too much
    % noise. Or we could derive the formula for acceleration ourself (or
    % maybe the Goldstein ref has it), that would be a good exercise.
    
    % see also http://www.u.arizona.edu/~pen/ame553/lessons.html, lesson 10
    % (see also the textbook at http://www.u.arizona.edu/~pen/ame553/) and 
    % http://www.rst.e-technik.tu-dortmund.de/cms/Medienpool/Downloads/Lehre/Vorlesungen/Robotics_Theory/Diff_Kinematics_4.pdf
else
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

if ETparams.data.qAlsoStoreComponentDerivs
    % also store velocities and acceleration in X and Y direction
    data.deg.velAz  = tempV(:,1) * ETparams.samplingFreq;
    data.deg.velEl  = tempV(:,2) * ETparams.samplingFreq;
    data.deg.accAz  = tempA(:,1) * ETparams.samplingFreq^2;
    data.deg.accEl  = tempA(:,2) * ETparams.samplingFreq^2;
end

if ETparams.data.qAlsoStoreandSmoothPixels
    % calculate derivative magnitudes
    data.pix.vel    = hypot(tempVpix(:,1), tempVpix(:,2)) * ETparams.samplingFreq;
    data.pix.acc    = hypot(tempApix(:,1), tempApix(:,2)) * ETparams.samplingFreq^2;
    
    % also store velocities and acceleration in X and Y direction
    if ETparams.data.qAlsoStoreComponentDerivs
        data.pix.velX   = tempVpix(:,1) * ETparams.samplingFreq;
        data.pix.velY   = tempVpix(:,2) * ETparams.samplingFreq;
        data.pix.accX   = tempApix(:,1) * ETparams.samplingFreq^2;
        data.pix.accY   = tempApix(:,2) * ETparams.samplingFreq^2;
    end
end





%%% helpers
function varargout = sgFilt(x,difforder,ntaps)
% wrapper for convenient syntax, fitting always using 2nd order polynomial
for p=1:length(difforder)
    varargout{p} = savitzkyGolayFilt(x,2,difforder(p),ntaps);
end