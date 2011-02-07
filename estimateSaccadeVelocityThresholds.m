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

% prepare algorithm parameters
minFixSamples       = ceil(ETparams.fixation.minDur/1000 * ETparams.samplingFreq);
centralFixSamples   = ceil(ETparams.saccade.minDur /6000 * ETparams.samplingFreq);

% assign initial thresholds
data.saccade.peakVelocityThreshold = ETparams.saccade.peakVelocityThreshold;
previousPeakDetectionThreshold = inf;

while previousPeakDetectionThreshold - data.saccade.peakVelocityThreshold > 1
    
    previousPeakDetectionThreshold = data.saccade.peakVelocityThreshold;
    
    % Find parts where the velocity is below the threshold, possible
    % fixation time (though this is still crude)
    qBelowThresh = data.deg.vel < data.saccade.peakVelocityThreshold;
    
    if nargin==2 || qusecentralsample
        % this is not just done to speed up iteration, we need to cut off
        % the edges of the testing intervals and only use parts of the data
        % that are likely to belong to fixations or the iteration will not
        % converge to a lower threshold. So always use this code path (just
        % call this function with 2 arguments), unless you want to see it
        % for yourself.
        
        % get bounds of these detected peaks
        [threshon,threshoff] = bool2bounds(qBelowThresh);
        
        % throw out intervals that are too short and therefore unlikely to
        % be fixations
        qLongEnough = threshoff-threshon >= minFixSamples;
        threshon = threshon (qLongEnough);
        threshoff= threshoff(qLongEnough);
        
        % shrink them as done in Nystrom's version, to make sure we don't
        % catch the data that is still during the saccade
        threshon = threshon +floor(centralFixSamples);
        threshoff= threshoff-ceil (centralFixSamples);
        
        % convert to data selection indices
        idx=bounds2ind(threshon,threshoff);
        
        % get mean and std of this data
        meanVel = nanmean(data.deg.vel(idx));
        stdVel  = nanstd (data.deg.vel(idx));
    else
        meanVel = nanmean(data.deg.vel(qBelowThresh));
        stdVel  = nanstd (data.deg.vel(qBelowThresh));
    end
    
    % calculate new thresholds
    data.saccade.peakVelocityThreshold  = meanVel + 6*stdVel;
    data.saccade.onsetVelocityTreshold  = meanVel + 3*stdVel;
end