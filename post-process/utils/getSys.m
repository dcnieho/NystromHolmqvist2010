function sys = getSys(order,dampening,qcontrol,qdisturb,qtracking)
% input the continuous-time state-space representation of the system,
% flexible and should be no issue after reading the document
assert(qcontrol||qdisturb,'no inputs requested')

if xor(qcontrol,qdisturb)
    % only one input, either control only or disturbance only
    assert(isscalar(order),'making system for one input only, but Yc orders for multiple inputs defined');
    
    [phi,gamma,cee,dee] = getMat(order,dampening);
    
elseif qtracking
    % tracking control setup: 2 inputs, 2 outputs
    % first  input is routed to first  output, and
    % second input is routed to second output
    assert(length(order)==2, '2 inputs so need two orders')
    % get system w.r.t. first input
    [phi,gamma,cee,dee] = getMat(order(1),dampening);
    % get system w.r.t. second input
    [sphi,sgamma,scee,sdee] = getMat(order(2),dampening);
    if order(2)==0 || order(1)==0
        gamma   = [gamma sgamma];
        cee     = [cee; scee];
    end
    if order(2)>0
        if order(1)==0
            % we need to use the phi matrix of the second input
            phi = sphi;
        else
            % more elaborate matrix structure needed to keep non-zeroth-order
            % disturbance input separate from the control input
            phi(3:4,3:4) = sphi;
            gamma(3:4,2) = sgamma;
            cee(2,3:4) = scee;
        end
    end
    dee     = [dee 0; 0 sdee];
    
else
    % compensatory control setup: 2 inputs, one output
    assert(length(order)==2, '2 inputs so need two orders')
    assert(order(2)<=order(1),'can''t do higher disturbance order than control order') % phi would be wrong with the current code... not sure if it could be fixed
    % get system w.r.t. first input
    [phi,gamma,cee,dee] = getMat(order(1),dampening);
    % get system w.r.t. second input
    [~,sgamma,~,sdee] = getMat(order(2),dampening);
    gamma   = [gamma sgamma];
    dee     = [dee sdee];
end

sys = ss(phi,gamma,cee,dee);



%%% helpers
function [phi,gamma,cee,dee] = getMat(order,dampening)

% we always have two states, even if 0th order system, just like in the
% document. Just for convenience that things are easy to compare
switch order
    case 0
        phi     = [0 0; 0 0];
        gamma   = [0; 0];
        cee     = [0 0];
        dee     = 1;
    case 1
        phi     = [0 0; 0 0];
        gamma   = [0; 1];
        cee     = [0 1];
        dee     = 0;
    case 2
        phi     = [-dampening 0; 1 0];
        gamma   = [1; 0];
        cee     = [0 1];
        dee     = 0;
    otherwise
        error('order %d not defined, but should be straightforward to do it yourself');
end