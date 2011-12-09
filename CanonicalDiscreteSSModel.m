function [out,state] = CanonicalDiscreteSSModel(sys,in)

[phi,gamma,cee,dee] = ssdata(sys);

if size(in,1) > size(in,2)
    in = in.';
end

datalen = length(in);
nstate = size(phi,1);
nout = size(cee,1);

% run through state transition matrices
state = zeros(nstate, datalen);
out   = zeros(nout  , datalen);

for p=1:length(in)
    state(:,p+1) = phi*state(:,p) + gamma*in(:,p);
    out(p)       = cee*state(:,p) + dee  *in(:,p);
end