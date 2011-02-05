function [data,qnoise] = detectAndRemoveNoise(data,ETparams)
% Detects and removes un-physiological movement (which derives from noise
% and blinks)

V_threshold = median(data.vel)*2;

% Detect possible blinks and noise (where XY-coords are 0 or if the eyes move too fast)
% TODO make check more general, checking if outside screen, current doesn't
% work for our EyeLink, which doesn't give 0 (though already nan in that
% case)
% do not have to process things that are already NaN
qnoise = (data.X <= 0 & data.Y <= 0) |...
        data.vel > ETparams.blinkVelocityThreshold |...
        abs(data.acc) > ETparams.blinkAccThreshold;

% Label blinks or noise
[noiseon,noiseoff] = findContiguousRegions(qnoise);

% Process one blink or noise period at the time, refine the bounds
% TODO: does this sometimes delete way too much data? (if go below
% threshold only near next noise or so.. dunno.. seems unlikely)
for p = 1:length(noiseon)

    % Go back in time to see where the blink (noise) started
    sEventIdx   = find(data.vel(noiseon(p):-1:1) <= V_threshold, 1);
    if isempty(sEventIdx), continue, end
    noiseon(p)  = noiseon(p) - sEventIdx + 1;
    
    % Go forward in time to see where the blink (noise) started    
    eEventIdx   = find(data.vel(noiseoff(p):end) <= V_threshold, 1);
    if isempty(eEventIdx), continue, end    
    noiseoff(p) = noiseoff(p) + eEventIdx - 1;
end

% create boolean matrix given refined noise bounds
qnoise = bounds2bool(noiseon,noiseoff,length(data.vel));
if sum(qnoise)/length(data.vel) > 0.20
    disp('Warning: This trial contains > 20 % noise+blinks samples')
    data.qNoiseTrial = true;
else
    data.qNoiseTrial = false;
end
% remove data that is due to noise
data.X(qnoise) = nan;
data.Y(qnoise) = nan;
data.vel(qnoise) = nan;
data.acc(qnoise) = nan;

% store noise bounds
data.noise = [noiseon noiseoff];