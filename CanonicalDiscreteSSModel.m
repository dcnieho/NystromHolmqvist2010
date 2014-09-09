function [out,state] = CanonicalDiscreteSSModel(sys,in)

if exist('lsim.m','file')==2
    Ts = getTs(sys);
    time = [0:Ts:(length(in)-1)*Ts];
    [out,~,state] = lsim(sys,in,time);
else
    [phi,gamma,cee,dee] = ssdata(sys);
    
    datalen = length(in);
    nstate = size(phi,1);
    nout = size(cee,1);
    
    % run through state transition matrices
    state = zeros(nstate, datalen);
    out   = zeros(nout  , datalen);
    
    for p=1:length(in)
        state(:,p+1) = phi*state(:,p) + gamma*in(p);
        out(p)       = cee*state(:,p) + dee  *in(p);
    end
    
    % prepare output
    state = state(:,1:end-1).';
    out   = out.';
end