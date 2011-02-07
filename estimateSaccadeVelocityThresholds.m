function data = estimateSaccadeVelocityThresholds(data,ETparams,qusecentralsample)
% iteratively establishes the best velocity threshold for saccade detection
% for this trial.
% Given a threshold T, we take all data of the trial where the eye velocity
% as less than T and calculate its mean and std. We then establish a new
% threshold at V = mean + 6*std. We repeat this until the new threshold is
% less than 1°/s lower than the old threshold.
% The thus established velocity threshold is clos to optimal in that it
% allows us to detect as many saccades as possible given the noise in the
% data but does not detect too many due to possibly high noise in the data.

% !! Call this function with two parameters, unless you know what you are
% doing.

% assign initial thresholds
data.peakDetectionThreshold = ETparams.peakDetectionThreshold;
previousPeakDetectionThreshold = inf;

while previousPeakDetectionThreshold - data.peakDetectionThreshold > 1
    
    previousPeakDetectionThreshold = data.peakDetectionThreshold;
    
    % Find parts where the velocity is below the threshold, possible
    % fixation time (though this is still crude)
    qBelowThresh = data.vel < data.peakDetectionThreshold;
    
    if nargin==2 || qusecentralsample
        % this is not just done to speed up iteration, we need to cut off
        % the edges of the testing intervals and only use parts of the data
        % that are likely to belong to fixations or the iteration will not
        % converge to a lower threshold. So always use this code path (just
        % call this function with 2 arguments), unless you want to see it
        % for yourself.
        
        % get bounds of these detected peaks
        [threshon,threshoff] = findContiguousRegions(qBelowThresh);
        
        % throw out intervals that are too short and therefore unlikely to
        % be fixations
        qLongEnough = (threshoff-threshon)./ETparams.samplingFreq >= ETparams.minFixDur;
        threshon = threshon (qLongEnough);
        threshoff= threshoff(qLongEnough);
        
        % shrink them as done in Nystrom's version, to make sure we don't
        % catch the data that is still during the saccade
        centralFixSamples = ETparams.minFixDur*ETparams.samplingFreq/6;
        threshon = threshon +floor(centralFixSamples);
        threshoff= threshoff-ceil (centralFixSamples);
        
        % convert to data selection indices
        idx=bounds2ind(threshon,threshoff);
        
        % get mean and std of this data
        meanVel = nanmean(data.vel(idx));
        stdVel  = nanstd (data.vel(idx));
    else
        meanVel = nanmean(data.vel(qBelowThresh));
        stdVel  = nanstd (data.vel(qBelowThresh));
    end
    
    % calculate new thresholds
    data.peakDetectionThreshold  = meanVel + 6*stdVel;
    data.saccadeVelocityTreshold = meanVel + 3*stdVel;
end