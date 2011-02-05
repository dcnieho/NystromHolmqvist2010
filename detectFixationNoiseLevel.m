function detectFixationNoiseLevel(i,j)
global ETparams
 
possibleFixationIdx = ~ETparams.data(i,j).InitialVelPeakIdx;
fixLabeled = bwlabel(possibleFixationIdx);

% Process one inter-peak-saccadic periods (called fixations below,
% although they are not identified as fixations yet). 
fixNoise = [];
for k = 1:max(fixLabeled)

    % The samples related to the current fixation
    fixIdx = find(fixLabeled == k);
    
    % Check that the fixation duration exceeds the minimum duration criteria. 
    if length(fixIdx)/ETparams.samplingFreq < ETparams.minFixDur
        continue    
    end
    
    % Extract the samples from the center of the fixation
    centralFixSamples = ETparams.minFixDur*ETparams.samplingFreq/6;
    fNoise = ETparams.data(i,j).vel(floor(fixIdx(1)+centralFixSamples):ceil(fixIdx(end)-centralFixSamples));
    fixNoise = [fixNoise fNoise];
end
% size(fixNoise)
ETparams.data(i,j).avgNoise = nanmean(fixNoise);
ETparams.data(i,j).stdNoise = nanstd(fixNoise);

% Base the peak velocity threshold on the noise level
ETparams.data(i,j).peakDetectionThreshold =  ETparams.data(i,j).avgNoise + 6*ETparams.data(i,j).stdNoise;
ETparams.data(i,j).saccadeVelocityTreshold = ETparams.data(i,j).avgNoise + 3*ETparams.data(i,j).stdNoise;
ETparams.data(i,j).velPeakIdx  = ETparams.data(i,j).vel > ETparams.data(i,j).peakDetectionThreshold;

