function data = processFixations(data,ETparams)
% Collects information about the fixations

%%% timing
% start time (in milliseconds)
data.fixation.start    = data.fixation.on / ETparams.samplingFreq * 1000;

% duration (in milliseconds)
data.fixation.duration = (data.fixation.off-data.fixation.on+1) / ETparams.samplingFreq * 1000;

%%% drift during fixation (defined as position at end minus position at
%%% begin)
% drift amplitude & direction
[data.fixation.driftAmplitude, data.fixation.driftDirection] = ...
    calcAmplitudeFick(...
        data.deg.Azi(data.fixation.on ), data.deg.Ele(data.fixation.on ),...
        data.deg.Azi(data.fixation.off), data.deg.Ele(data.fixation.off)...
    );

       
%%% and now some that are best done in a for-loop
data.fixation.meanX_deg         = zeros(size(data.fixation.on));
data.fixation.meanY_deg         = zeros(size(data.fixation.on));
data.fixation.meanX_pix         = zeros(size(data.fixation.on));
data.fixation.meanY_pix         = zeros(size(data.fixation.on));
data.fixation.meanVelocity      = zeros(size(data.fixation.on));
data.fixation.peakVelocity      = zeros(size(data.fixation.on));
data.fixation.meanAcceleration  = zeros(size(data.fixation.on));
data.fixation.peakAcceleration  = zeros(size(data.fixation.on));
data.fixation.RMSS2S            = zeros(size(data.fixation.on));
data.fixation.STDx              = zeros(size(data.fixation.on));
data.fixation.STDy              = zeros(size(data.fixation.on));
data.fixation.BCEA              = zeros(size(data.fixation.on));
for p=1:length(data.fixation.on)
    idxs = data.fixation.on(p) : data.fixation.off(p);
    
    % average eye position
    data.fixation.meanX_deg(p)          = nanmean(data.deg.Azi(idxs));
    data.fixation.meanY_deg(p)          = nanmean(data.deg.Ele(idxs));
    % and convert it to pixels
    [data.fixation.meanX_pix(p),data.fixation.meanY_pix(p)] = ...
        fick2pix(data.fixation.meanX_deg(p), data.fixation.meanY_deg(p),ETparams);
    
    % mean and peak velocity and acceleration
    data.fixation.meanVelocity(p)       = nanmean(data.deg.vel(idxs));
    data.fixation.peakVelocity(p)       = max    (data.deg.vel(idxs));
    data.fixation.meanAcceleration(p)   = nanmean(data.deg.acc(idxs));
    data.fixation.peakAcceleration(p)   = max    (data.deg.acc(idxs));
    
    % fixation instability
    % calculate RMS S2S
    % since its done with diff, don't just exclude missing and treat
    % resulting as one continuous vector. replace missing with nan first,
    % use left-over values
    xdif = subsasgn(data.deg.Azi(idxs),substruct('()',{isnan(data.deg.vel(idxs))}),nan);
    ydif = subsasgn(data.deg.Ele(idxs),substruct('()',{isnan(data.deg.vel(idxs))}),nan);
    xdif = diff(xdif).^2; xdif(isnan(xdif)) = [];
    ydif = diff(ydif).^2; ydif(isnan(ydif)) = [];
    data.fixation.RMSS2S(p) = sqrt(mean(xdif + ydif)); % Sample 2 sample displacement RMS
    
    % calculate STD
    xposf= data.deg.Azi(idxs);
    yposf= data.deg.Ele(idxs);
    data.fixation.STDx(p) = std(xposf(~qMiss));
    data.fixation.STDy(p) = std(yposf(~qMiss));
    
    % calculate BCEA (Crossland and Rubin 2002 Optometry and Vision Science)
    qMiss= isnan(data.deg.vel(idxs));
    xx   = corrcoef(xposf(~qMiss),yposf(~qMiss));
    rho  = xx(1,2);
    P    = 0.68; % cumulative probability of area under the multivariate normal
    k    = log(1/(1-P));
    
    data.fixation.BCEA(p) = 2*k*pi*data.fixation.STDx(p)*data.fixation.STDy(p)*sqrt(1-rho.^2);
end
