function data = processFixations(data,ETparams)
% Collects information about the fixations

%%% timing
% start time (in milliseconds)
data.fixation.start    = data.fixation.on / ETparams.samplingFreq;

% duration (in milliseconds)
data.fixation.duration = (data.fixation.off-data.fixation.on) / ETparams.samplingFreq;

%%% drift during fixation (defined as position at end minus position at
%%% begin)
% drift amplitude & direction
[data.fixation.driftAmplitude, data.fixation.driftDirection] = ...
    calcAmplitudeFick(...
        data.deg.X(data.fixation.on ), data.deg.Y(data.fixation.on ),...
        data.deg.X(data.fixation.off), data.deg.Y(data.fixation.off)...
    );

       
%%% and now some that are best done in a for-loop
data.fixation.meanX             = zeros(size(data.fixation.on));
data.fixation.meanY             = zeros(size(data.fixation.on));
data.fixation.meanVelocity      = zeros(size(data.fixation.on));
data.fixation.peakVelocity      = zeros(size(data.fixation.on));
data.fixation.meanAcceleration  = zeros(size(data.fixation.on));
data.fixation.peakAcceleration  = zeros(size(data.fixation.on));
for p=1:length(data.fixation.on)
    idxs = data.fixation.on(p) : data.fixation.off(p);
    
    % average eye position
    data.fixation.meanX(p)              = mean(data.deg.X(idxs));
    data.fixation.meanY(p)              = mean(data.deg.Y(idxs));
    
    % mean and peak velocity and acceleration
    data.fixation.meanVelocity(p)       = mean(data.deg.vel(idxs));
    data.fixation.peakVelocity(p)       = max (data.deg.vel(idxs));
    data.fixation.meanAcceleration(p)   = mean(data.deg.acc(idxs));
    data.fixation.peakAcceleration(p)   = max (data.deg.acc(idxs));
end