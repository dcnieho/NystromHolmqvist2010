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

% find bounds of blinks or noise as detected above
[noiseon,noiseoff] = findContiguousRegions(qnoise);

% find bounds of data above threshold
[threshon,threshoff] = findContiguousRegions(data.vel > V_threshold);

% Process one blink or noise period at the time, refine the bounds
% We refine using the velocity threshold crosses computed on top. As this
% threshold is lower than the blink velocity threshold, on- and offset of
% this (very low!) median velocity threshold should enclose blink velocity
% thresholds, but not necessarily for any of the other sources of noise.
% Only if the noise bounds are enclosed by the median velocity bounds,
% replace. This also removes the possibility that we latch on to a wrong
% bound much further away in the data, as was possible in the old version
% (although not likely to generate much error).
for p = 1:length(noiseon)
    qenclosed = threshon<=noiseon(p) & threshoff>=noiseoff(p);
    
    if any(qenclosed)
        noiseon(p)  = threshon (qenclosed);
        noiseoff(p) = threshoff(qenclosed);
    end
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