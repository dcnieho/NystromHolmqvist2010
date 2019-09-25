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
    
% diff pupil size. always use Savitzky-Golay as some smoothing is needed to
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
% Calculate eye velocity and acceleration straightforwardly by applying
% Pythagoras' theorem. This gives us no information about the
% instantaneous axis of the eye rotation, but eye velocity is
% calculated correctly. Apply scale for velocity, as a 10° azimuth
% rotation at 0° elevation does not cover same distance as it does at
% 45° elevation: sqrt(theta_dot^2*cos^2 phi + phi_dot^2)
if 0
    % this is as per Warren, R. (1982). Optical transformations during
    % movement: review of the optical concomitants of egospeed. Report
    % AFOSR-TR-82-1028, Columbus: Ohio State University, Department of
    % Psychology, Aviation Psychology Laboratory.
    data.deg.vel    = hypot(tempV(:,1).*cosd(data.deg.Ele), tempV(:,2));
    % No idea what scaling to apply for acceleration, figure that out
    % some time...
    data.deg.acc    = hypot(tempA(:,1), tempA(:,2));
else
    % use method as per Singh, T., Perry, C., & Herter, T. (2015). A
    % geometric method for computing ocular kinematics and classifying
    % gaze events using monocular remote eye tracking in a robotic
    % environment. Journal of Neuroengineering and Rehabilitation.
    % This uses a different spherical coordinate system, but that
    % doesn't matter as we the final outcome of angular velocity
    % computations is the same, and they also provide a computation for
    % acceleration!
    % 1. prep input
    pixPerMeter = ETparams.screen.resolution ./ ETparams.screen.size;
    px          = data.pix.X./pixPerMeter(1);
    py          = data.pix.Y./pixPerMeter(2);
    pz          = repmat(ETparams.screen.viewingDist,size(px));
    hpxpy       = hypot(px,py);
    rho         = hypot(hpxpy,pz);
    theta       = atan2(py,px);
    phi         = acos(pz./rho);
    % 2. calculate derivatives
    pxd         = -sgFilt(px,1,ntaps) * ETparams.samplingFreq;
    pyd         = -sgFilt(py,1,ntaps) * ETparams.samplingFreq;
    pzd         = 0;    % derivative of constant is 0, so just hardcode that
    % 3. calculate component velocities as per eqs. 4a and 4b in Singh
    % et al. (2015)
    theta_d     = (pxd.*py-px.*pyd)./(px.^2+py.^2);
    phi_d       = pz.*(px.*pxd+py.*pyd)./(rho.^2.*sqrt(px.^2+py.^2));       % simplified from equation in paper because pzd==0
    rho_d       = (px.*pxd+py.*pyd+pz.*pzd)./rho;
    % 4. calculate angular velocity
    data.deg.vel= hypot(theta_d.*sin(phi),phi_d)*180/pi;
    % 5. calculate angular acceleration
    theta_dd    = conv2(theta_d,[1 0 -1].','valid')/2 * ETparams.samplingFreq; theta_dd = theta_dd([1 1:end end]);
    phi_dd      = conv2(  phi_d,[1 0 -1].','valid')/2 * ETparams.samplingFreq;   phi_dd =   phi_dd([1 1:end end]);
    thetaTerm   = 2.*rho_d.*theta_d.*sin(phi)./rho+theta_dd.*sin(phi)-2.*theta_d.*phi_d.*cos(phi);
    phiTerm     = 2.*rho_d.*phi_d./rho+phi_d.^2.*sin(phi).*cos(phi)+phi_dd;
    data.deg.acc= hypot(thetaTerm,phiTerm)*180/pi;
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
