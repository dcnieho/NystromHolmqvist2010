function data = filterDataAndCalcDerivative(data,ETparams)
% By default, we'll use the Savitzky-Golay filter for smoothing and
% differentiation (simply numerical difference calculation is also
% possible), which basically performs a polynomial least square fit in a
% moving window. The eye position data is smoothed by extracting the zeroth
% order coefficient of this polynomial. Eye velocity and acceleration and
% obtained by returning the higher order coefficients of the fitted
% polynomial. Note that for kth order derivatives (k>1) we'd need to
% multiply the ouput by k! to obtain the derivative coefficient (see Taylor
% series expansion). This is already taken into account internally in the
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
% it. I am not sure about the extent to which this increases the number of
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
% In general, we'll need to assess how problematic this really is. When it
% comes to the filter, we need to simulate its performance (how much
% widening, how much undershoot is there really?). Generate a function much
% like the eye movement, but without any noise (Gaussian velocity profile,
% very low std for sharp peak? Or more directly cut a part from some
% eyemovement data and smooth it by hand) and see what happens to it when
% filtering by polynomial fit... Also, remember the Luo 2005 paper on
% properties of Savitzky-Golay differentiation filters. By the way, i
% indeed think my logic above was correct: the zeroth order moment of the
% first derivative of a time series should correspond/be relatable to the
% first order moment of the time series itself. It is thus a problem that
% the width of the velocity peak (3rd moment) is not present by a second
% order polynomial fit.
%
% Also look at:
% filtering / digital differentiation for the eye data, look at Pavel
% Holoborodko's stuff: http://www.holoborodko.com/pavel/numerical-methods/.
% Not published I think, but it might be worth it, and it seems possible to
% contact the guy, though he might not respond quickly...
% http://nl.wikipedia.org/wiki/Savitsky-Golay_filter


% calculate component velocities and accelerations
%--------------------------------------------------------------------------
if ETparams.data.qNumericallyDifferentiate
    if ETparams.data.qNumericallyDifferentiate==1
        % instead of diff we use conv2(data,[1 0 -1].','valid')/2 because it is centered correctly
        % NB diff() is equivalent to conv2(data,[1 -1].','valid')
        tempV   = conv2([data.deg.Azi data.deg.Ele],[1 0 -1].','valid')/2 * ETparams.samplingFreq;
        tempV   = tempV([1 1:end end],:);    % deal with end effects
        % diff always leads to one extra sample of missing at both sides,
        % compensate for this in the same way as we deal with other data
        % end effects
        [on,off] = bool2bounds(isnan(tempV(:,1)));
        for p=1:length(on)
            if on(p)>1
                tempV(on (p),:) = tempV(on (p)-1,:);
            end
            if off(p)<size(tempV,1)
                tempV(off(p),:) = tempV(off(p)+1,:);
            end
        end
        
        tempA   = conv2(tempV,[1 0 -1].','valid')/2 * ETparams.samplingFreq^2;
        tempA   = tempA([1 1:end end],:);    % deal with end effects
        [on,off] = bool2bounds(isnan(tempA(:,1)));  % deal with end effects due to missing
        for p=1:length(on)
            if on(p)>1
                tempA(on (p),:) = tempA(on (p)-1,:);
            end
            if off(p)<size(tempA,1)
                tempA(off(p),:) = tempA(off(p)+1,:);
            end
        end
        
        % also calculate derivatives for eye position in pixels
        if ETparams.data.qAlsoStoreandDiffPixels
            tempVpix   = conv2([data.pix.X data.pix.Y],[1 0 -1].','valid')/2 * ETparams.samplingFreq;
            tempVpix   = tempVpix([1 1:end end],:);    % deal with end effects
            
            tempApix   = conv2(tempVpix,[1 0 -1].','valid')/2 * ETparams.samplingFreq^2;
            tempApix   = tempApix([1 1:end end],:);    % deal with end effects
        end
    elseif ETparams.data.qNumericallyDifferentiate==2
        % simple diff
        tempV   = diff([data.deg.Azi data.deg.Ele],1,1);
        tempA   = diff([data.deg.Azi data.deg.Ele],2,1);
        
        % make same length as position trace by repeating last sample
        tempV   = tempV([1:end end    ],:) * ETparams.samplingFreq;
        tempA   = tempA([1:end end end],:) * ETparams.samplingFreq^2;
        
        % also calculate derivatives for eye position in pixels
        if ETparams.data.qAlsoStoreandDiffPixels
            tempVpix   = diff([data.pix.X data.pix.Y],1,1);
            tempApix   = diff([data.pix.X data.pix.Y],2,1);
            
            % make same length as position trace by repeating last sample
            tempVpix   = tempVpix([1:end end    ],:) * ETparams.samplingFreq;
            tempApix   = tempApix([1:end end end],:) * ETparams.samplingFreq^2;
        end
    end
else
    % prepare parameters
    %--------------------------------------------------------------------------
    % span of filter, use minimum length of saccade. Its very important to not
    % make the filter window wider than the narrowest feature we are interested
    % in, or we'll smooth out those features.
    window  = ceil(ETparams.data.filterWindow/1000*ETparams.samplingFreq);
    % number of filter taps
    ntaps   = 2*ceil(window)-1;
    
    % calculate derivatives
    [tempV,tempA] = sgFilt([data.deg.Azi data.deg.Ele],[1 2],ntaps);
    tempV = -tempV * ETparams.samplingFreq;                 % not sure why, but the Savitzky-Golay filter gives me the wrong sign for the component velocities
    tempA =  tempA * ETparams.samplingFreq^2;               % note that no need to multiply by factorial(2) as filter coefficients used already include this scaling
    
    % also calculate derivatives for eye position in pixels
    if ETparams.data.qAlsoStoreandDiffPixels
        [tempVpix,tempApix] = sgFilt([data.pix.X data.pix.Y],[1 2],ntaps);
        tempVpix = -tempVpix * ETparams.samplingFreq;       % not sure why, but the Savitzky-Golay filter gives me the wrong sign for the component velocities
        tempApix =  tempApix * ETparams.samplingFreq^2;     % note that no need to multiply by factorial(2) as filter coefficients used already include this scaling
    end
    
end
    
% diff pupil size. always use savitzky golay as some smoothing is needed to
% use this change of pupil size in a meaningful way for blink detection
if isfield(data,'pupil') && ~isempty(data.pupil.size)
    if ETparams.data.qNumericallyDifferentiate==1 && 0
        data.pupil.dsize = conv2(data.pupil.size,[1 0 -1].','valid')/2 * ETparams.samplingFreq;
        data.pupil.dsize = [data.pupil.dsize(1,:); data.pupil.dsize; data.pupil.dsize(end,:)];    % deal with end effects
    elseif ETparams.data.qNumericallyDifferentiate==2 && 0
        data.pupil.dsize = diff(data.pupil.size,1,1) * ETparams.samplingFreq;
        data.pupil.dsize = data.pupil.dsize([1:end end],:);    % deal with end effects
    else
        window  = ceil(ETparams.data.filterWindow/1000*ETparams.samplingFreq);
        ntaps   = 2*ceil(window)-1;
        data.pupil.dsize = -sgFilt(data.pupil.size,1,ntaps) * ETparams.samplingFreq;
    end
end

% compute eye velocities, and accelerations
%--------------------------------------------------------------------------
if ETparams.data.qPreciseCalcDeriv
    % TODO: implement, acceleration still TODO
    % also TODO: check signs! For Haslwanter, positive corresponds to
    % leftward azimuth, downward elevation and clockwise torsion. As noted
    % in prepareData and below, our axes do not form a consistent
    % right-handed system and should thus be calculating the wrong axis
    % here...
    % revise our coordinate system first when needed!
    error('TODO')
    
    % Now calculate eye velocity and acceleration precisely
    % See Equation 30 in Haslwanter T (1995) Mathematics of 3-dimensional
    % eye rotations, Vision Res 35, 1727-1739. (look up that original
    % Goldstein reference as well).
    % This allows us to use the analytical differentiation of the
    % polynomial fitted to eye azimuth and elevation (so-called coordinate
    % velocities).
    % The eye velocity vector is commonly called omega, where the magnitude
    % of omega is the angular velocity and its axis indicates the
    % instantaneous axis of the eye rotation. Omega, a vector in the
    % head-fixed reference frame, is ordered Z-X-Y where a vector the Z
    % direction points in the primary position, an X vector points leftward
    % and an Y vector upward. A rotation along Z is roll/torsion, along X
    % is pitch/elevation and along Y is yaw/azimuth (see Haslwanter 1995,
    % Eq. 1 and Figure 1). As noted in prepareData, our current head
    % reference axes do not conform nor even form a consistently
    % righthanded system. When interested in using this code, revise our
    % axes!
    
    % see also http://www.u.arizona.edu/~pen/ame553/lessons.html, lesson 10
    % (see also the textbook at http://www.u.arizona.edu/~pen/ame553/)
    % -> has acceleration as a quaternion! Textbook (appendix) has
    % equivalence formulas for velocity as quaternion and other schemes if
    % i remember correctly.
    % also: 
    % http://www.rst.e-technik.tu-dortmund.de/cms/Medienpool/Downloads/Lehre/Vorlesungen/Robotics_Theory/Diff_Kinematics_4.pdf
    
    % formula 30 in Haslwanter is for 3D rotation, but we have no knowledge
    % of the torsional component, so we'll have to do with that set to 0
    if 0
        % if we had 3D data, and tempV(:,3) is the torsional velocity:
        data.deg.omega = [...
            tempV(:,3).*cosd(data.deg.Azi).*cosd(data.deg.Ele) - tempV(:,2).*sind(data.deg.Azi)                     , ...
            tempV(:,2).*cosd(data.deg.Azi)                     + tempV(:,3).*sind(data.deg.Azi).*cosd(data.deg.ele) , ...
            tempV(:,1)                                         - tempV(:,3).*sind(data.deg.Ele)                       ...
            ];
    else
        % for lack of information about torsion, we set those parts to 0
        % this means the instantaneous axis does not take changes in
        % torsion into account and might therefore be systematically off.
        % NOTE, that's a crap idea. below in else block has correct way to
        % do velocity in Fick system
        data.deg.omega = [...
                                                               - tempV(:,2).*sind(data.deg.Azi)                     , ...
            tempV(:,2).*cosd(data.deg.Azi)                                                                          , ...
            tempV(:,1)                                                                                                ...
            ];
    end
    data.deg.vel    = sqrt(sum(data.deg.omega.^2,2));
else
    % Calculate eye velocity and acceleration straightforwardly by applying
    % Pythagoras' theorem. This gives us no information about the
    % instantaneous axis of the eye rotation, but eye velocity is
    % calculated correctly. Apply scale for velocity, as a 10° azimuth
    % rotation at 0° elevation does not cover same distance as it does at
    % 45° elevation: sqrt(theta_dot^2*cos^2 phi + phi_dot^2)
    % No idea what scaling to apply for acceleration, figure that out some
    % time...
    data.deg.vel    = hypot(tempV(:,1).*cosd(data.deg.Ele), tempV(:,2));
    data.deg.acc    = hypot(tempA(:,1), tempA(:,2));
end

if ETparams.data.qAlsoStoreComponentDerivs
    % also store velocities and acceleration in X and Y direction
    data.deg.velAzi = tempV(:,1);
    data.deg.velEle = tempV(:,2);
    data.deg.accAzi = tempA(:,1);
    data.deg.accEle = tempA(:,2);
end

if ETparams.data.qAlsoStoreandDiffPixels
    % calculate derivative magnitudes
    data.pix.vel    = hypot(tempVpix(:,1), tempVpix(:,2));
    data.pix.acc    = hypot(tempApix(:,1), tempApix(:,2));
    
    % also store velocities and acceleration in X and Y direction
    if ETparams.data.qAlsoStoreComponentDerivs
        data.pix.velX   = tempVpix(:,1);
        data.pix.velY   = tempVpix(:,2);
        data.pix.accX   = tempApix(:,1);
        data.pix.accY   = tempApix(:,2);
    end
end





%%% helpers
function varargout = sgFilt(x,difforder,ntaps)
% wrapper for convenient syntax, fitting always using 2nd order polynomial

varargout = cell(1,length(difforder));
for p=1:length(difforder)
    varargout{p} = savitzkyGolayFilt(x,2,difforder(p),ntaps);
end
