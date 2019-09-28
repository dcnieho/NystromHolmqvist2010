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
if ETparams.data.numericallyDifferentiate
    if ETparams.data.numericallyDifferentiate==1
        % instead of diff we use conv2(data,[1 0 -1].','valid')/2 because it is centered correctly
        % NB diff() is equivalent to conv2(data,[1 -1].','valid')
        tempV   = conv2([data.deg.Azi data.deg.Ele data.deg.dist],[1 0 -1].','valid')/2 * ETparams.samplingFreq;
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
        if ETparams.data.alsoStoreandDiffPixels
            tempVpix   = conv2([data.pix.X data.pix.Y],[1 0 -1].','valid')/2 * ETparams.samplingFreq;
            tempVpix   = tempVpix([1 1:end end],:);    % deal with end effects
            
            tempApix   = conv2(tempVpix,[1 0 -1].','valid')/2 * ETparams.samplingFreq^2;
            tempApix   = tempApix([1 1:end end],:);    % deal with end effects
        end
    elseif ETparams.data.numericallyDifferentiate==2
        % simple diff
        tempV   = diff([data.deg.Azi data.deg.Ele data.deg.dist],1,1);
        tempA   = diff([data.deg.Azi data.deg.Ele data.deg.dist],2,1);
        
        % make same length as position trace by repeating last sample
        tempV   = tempV([1:end end    ],:) * ETparams.samplingFreq;
        tempA   = tempA([1:end end end],:) * ETparams.samplingFreq^2;
        
        % also calculate derivatives for eye position in pixels
        if ETparams.data.alsoStoreandDiffPixels
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
    [tempV,tempA] = sgFilt([data.deg.Azi data.deg.Ele data.deg.dist],[1 2],ntaps);
    tempV = -tempV * ETparams.samplingFreq;                 % not sure why, but the Savitzky-Golay filter gives me the wrong sign for the component velocities
    tempA =  tempA * ETparams.samplingFreq^2;               % note that no need to multiply by factorial(2) as filter coefficients used already include this scaling
    
    % also calculate derivatives for eye position in pixels
    if ETparams.data.alsoStoreandDiffPixels
        [tempVpix,tempApix] = sgFilt([data.pix.X data.pix.Y],[1 2],ntaps);
        tempVpix = -tempVpix * ETparams.samplingFreq;       % not sure why, but the Savitzky-Golay filter gives me the wrong sign for the component velocities
        tempApix =  tempApix * ETparams.samplingFreq^2;     % note that no need to multiply by factorial(2) as filter coefficients used already include this scaling
    end
    
end
    
% diff pupil size. always use Savitzky-Golay as some smoothing is needed to
% use this change of pupil size in a meaningful way for blink
% classification
if isfield(data,'pupil') && ~isempty(data.pupil.size)
    if ETparams.data.numericallyDifferentiate==1 && 0
        data.pupil.dsize = conv2(data.pupil.size,[1 0 -1].','valid')/2 * ETparams.samplingFreq;
        data.pupil.dsize = [data.pupil.dsize(1,:); data.pupil.dsize; data.pupil.dsize(end,:)];    % deal with end effects
    elseif ETparams.data.numericallyDifferentiate==2 && 0
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
% Calculate eye velocity and acceleration as length of their vectors.
% Derivation of the below formulas is in the subfunction deriveVelAcc at
% the bottom of this file

% Velocity: sqrt(theta_dot^2*cos^2 phi + phi_dot^2)
% applied scaling factor for velocity is intuitive: need to take into
% account that a 10° azimuth rotation at 0° elevation does not cover same
% distance as it does at 45° elevation
data.deg.vel    = hypot(tempV(:,1).*cosd(data.deg.Ele), tempV(:,2));
% Acceleration:
%     /  /                    2          2 Rd phd \ 2    /                                  2 Rd thd cos(ph) \ 2 \ 
% sqrt|  | cos(ph) sin(ph) thd  + phdd + -------- |   +  | thdd cos(ph) - phd thd sin(ph) + ---------------- |   | 
%     \  \                                   R    /      \                                          R        /   / 
thetaPart       = tempA(:,1).*cos(data.deg.Ele)-tempV(:,2).*tempV(:,1).*sin(data.deg.Ele)+(tempV(:,3).*tempV(:,1).*cos(data.deg.Ele).*2.0)./data.deg.dist;
phiPart         = tempA(:,2)+tempV(:,1).^2.*cos(data.deg.Ele).*sin(data.deg.Ele)+(tempV(:,3).*tempV(:,2).*2.0)./data.deg.dist;
data.deg.acc    = hypot(thetaPart,phiPart);

if ETparams.data.alsoStoreComponentDerivs
    % also store velocities and acceleration in X and Y direction
    data.deg.velAzi = tempV(:,1);
    data.deg.velEle = tempV(:,2);
    data.deg.velDist= tempV(:,3);
    data.deg.accAzi = tempA(:,1);
    data.deg.accEle = tempA(:,2);
    data.deg.accDist= tempA(:,3);
end

if ETparams.data.alsoStoreandDiffPixels
    % calculate derivative magnitudes
    data.pix.vel    = hypot(tempVpix(:,1), tempVpix(:,2));
    data.pix.acc    = hypot(tempApix(:,1), tempApix(:,2));
    
    % also store velocities and acceleration in X and Y direction
    if ETparams.data.alsoStoreComponentDerivs
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


function [velfunc, accfunc] = derivVelAcc
% this derives formulas for angular velocity and acceleration from
% spherical component derivatives.
% Requires symbolic math toolbox
% tested with MATLAB R2019a
syms R Rh Rd Rhd Rdd th thh thhd thdd ph phd phh phhd phdd
syms fR(t) fRh(t) fRhd(t) fth(t) fthd(t) fthh(t) fph(t) fphd(t) fphh(t)
% two variable names that shadow existing function names, need to create
% like this:
thd = sym('thd');
fRd = str2sym('fRd(t)');
if 1
    % Fick coordinates
    
    % 1. get transform Cartesian -> spherical
    % MATLAB's cart2sph models a Fick gimbal, although their reference Z axis
    % is our Y axis, their X axis is our Z axis and their Y axis is our X axis:
    % cart2sph: X Y Z
    % our Fick: Z X Y
    [z,x,y] = sph2cart(th,ph,R);
    sphRep  = [x y z]
    % 2. get unit vectors of this system
    rh_h    = diff(sphRep,R);               % this is already a unit vector
    th_h    = diff(sphRep,th)/R/cos(ph);    % normalize
    ph_h    = diff(sphRep,ph)/R;            % normalize
    % put together into transform matrix, turn into unit vecs
    Mf      = [rh_h;th_h;ph_h]
    uVecs   = sum(Mf,2);
    % output:
    % [ rho_hat ]                                                   [x_hat]
    % [theta_hat] =                      Mf                       * [y_hat]
    % [ phi_hat ]                                                   [z_hat]
    %               [ cos(ph)*sin(th), sin(ph),  cos(ph)*cos(th)]   [x_hat]
    %             = [         cos(th),       0,         -sin(th)] * [y_hat]
    %               [-sin(ph)*sin(th), cos(ph), -cos(th)*sin(ph)]   [z_hat]
    % this is an orthogonal matrix, because
    % simplify(inv(Mf),100)
    % equals
    % transpose(Mf)
    assert(all(all(logical(eval(simplify(inv(Mf),100)==transpose(Mf))))))
    % 3. get partial derivatives of unit vectors, will need them below
    uVecPartDerivs = jacobian(uVecs,[R, th, ph]);
    % 3.1 simplify by substituting in equations for unit vectors
    uVecPartDerivs = subs(uVecPartDerivs,uVecs.',[Rh thh phh]);
    % do some manually
    assert(logical(simplify(uVecPartDerivs(1,2)==uVecs(2)*cos(ph))))
    uVecPartDerivs(1,2) = thh*cos(ph);
    assert(logical(simplify(uVecPartDerivs(3,2)==-uVecs(2)*sin(ph))))
    uVecPartDerivs(3,2) = -thh*sin(ph);
    assert(logical(simplify(uVecPartDerivs(2,2)==-uVecs(1)*cos(ph) +uVecs(3)*sin(ph))))
    uVecPartDerivs(2,2) = -Rh*cos(ph)+phh*sin(ph)

else
    % standard (?), e.g. system as used on wikipedia etc
    % derivations follow (and match! yay!)
    % https://www.cpp.edu/~ajm/materials/delsph.pdf
    
    % 1. get transform Cartesian -> spherical
    x = R*sin(th)*cos(ph);
    y = R*sin(th)*sin(ph);
    z = R*cos(th);
    sphRep  = [x y z]
    % 2. get unit vectors of this system
    rh_h    = diff(sphRep,R);               % this is already a unit vector
    th_h    = diff(sphRep,th)/R;            % normalize
    ph_h    = diff(sphRep,ph)/R/sin(th);    % normalize
    % put together into transform matrix, turn into unit vecs
    Ms      = [rh_h;th_h;ph_h]
    uVecs   = sum(Ms,2);
    
    uVecPartDerivs = jacobian(uVecs,[R, th, ph]);
    uVecPartDerivs = subs(uVecPartDerivs,uVecs.',[Rh thh phh]);
    % do some manually
    assert(logical(simplify(uVecPartDerivs(1,3)==uVecs(3)*sin(th))))
    uVecPartDerivs(1,3) = phh*sin(th);
    assert(logical(simplify(uVecPartDerivs(2,3)==uVecs(3)*cos(th))))
    uVecPartDerivs(2,3) = phh*cos(th);
    assert(logical(simplify(uVecPartDerivs(3,3)==-uVecs(1)*sin(th) -uVecs(2)*cos(th))))
    uVecPartDerivs(3,3) = -Rh*sin(th) -thh*cos(th)
end

% get unit vector time derivatives
uVecDots = sum(repmat([Rd thd phd],3,1).*uVecPartDerivs,2)
% make into functional version for derivations below
fuVecDot1= subs(uVecDots(1),[phd phh thd thh ph th],[fphd fphh fthd fthh fph fth]);
% 4. now we can derive angular velocity:
pos     = fR(t)*fRh(t);                                                 % position of element in spherical coordinate system, rho*rho_hat
fvelvec = diff(pos,t);                                                  % diff position to get velocity
fvelvec = subs(fvelvec,[diff(fRh,t) diff(fR,t)],[fRhd fRd])             % substitute diff with the symbolic functions we have for those
velvec  = expand(subs(fvelvec,fRhd,uVecDots(1)))                        % now substitute derivative of Rh for Rhd
fvelvec = expand(subs(fvelvec,fRhd,fuVecDot1));                         % do same with functions version, need this for acceleration below
velvec  = children(velvec);                                             % split into components along unit vectors
velvec  = subs(subs(velvec/fR,fRh,0), [thh phh], [1 1]);                % remove unit vectors from formula (keep only tangential components), divide by R to remove it (can't explain exactly why)
pretty(norm(velvec))

% 5. accel
faccvec = diff(fvelvec,t);                                              % diff velocity to get acceleration
% replace all diffs with symbolic equivalents
faccvec = subs(faccvec,[diff(fRh,t) diff(fRd,t) diff(fthh,t) diff(fthd,t) diff(fR,t) diff(fphh,t), diff(fphd,t) diff(fth,t)],[Rhd Rdd thhd thdd Rd phhd phdd thd]);
% replace all functionals with symbolic equivalents
accvec  = subs(faccvec,[fR fRd fRh fph fphd fphh fth fthd fthh],[R Rd Rh ph phd phh th thd thh])
% we can substitute some of these symbols (differentials of unit vectors)
accvec  = subs(subs(subs(accvec,Rhd,sum(uVecDots(1,:))),thhd,sum(uVecDots(2,:))),phhd,sum(uVecDots(3,:)));
accvec  = collect(collect(collect(accvec,phh),thh),Rh);
% split into components along unit vectors
accvec  = children(accvec);
accvec  = [accvec(1) children(accvec(2))].';
% remove unit vectors from formula (keep only tangential components)
accvec  = subs(accvec,[Rh thh phh],[0 1 1]);
% need to help a bit with dividing R into it (can't explain exactly why need to divide R into it)
accvec(2)   = sum(children(accvec(2))/R);
accvec(3)   = sum(children(accvec(3))/R)
pretty(norm(accvec))


velfunc = matlabFunction(norm(velvec));
accfunc = matlabFunction(norm(accvec));