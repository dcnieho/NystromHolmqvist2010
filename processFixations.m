function data = processFixations(data,ETparams)
% Collects information about the fixations

%%% timing
% start time (in milliseconds)
data.fixation.start    = data.fixation.on / ETparams.samplingFreq;

% duration (in milliseconds)
data.fixation.duration = (data.fixation.off-data.fixation.on+1) / ETparams.samplingFreq;

%%% drift during fixation (defined as position at end minus position at
%%% begin)
% drift amplitude & direction
[data.fixation.driftAmplitude, data.fixation.driftDirection] = ...
    calcAmplitudeFick(...
        data.deg.X(data.fixation.on ), data.deg.Y(data.fixation.on ),...
        data.deg.X(data.fixation.off), data.deg.Y(data.fixation.off)...
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
for p=1:length(data.fixation.on)
    idxs = data.fixation.on(p) : data.fixation.off(p);
    
    % average eye position
    data.fixation.meanX_deg(p)          = nanmean(data.deg.X(idxs));
    data.fixation.meanY_deg(p)          = nanmean(data.deg.Y(idxs));
    % and convert it to degrees
    [data.fixation.meanX_pix(p),data.fixation.meanY_pix(p)] = ...
        fick2pix(data.fixation.meanX_deg(p), data.fixation.meanY_deg(p),ETparams);
    
    % mean and peak velocity and acceleration
    data.fixation.meanVelocity(p)       = nanmean(data.deg.vel(idxs));
    data.fixation.peakVelocity(p)       = max    (data.deg.vel(idxs));
    data.fixation.meanAcceleration(p)   = nanmean(data.deg.acc(idxs));
    data.fixation.peakAcceleration(p)   = max    (data.deg.acc(idxs));
end