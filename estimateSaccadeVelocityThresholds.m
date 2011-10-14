function data = estimateSaccadeVelocityThresholds(data,ETparams,qusecentralsample)
% iteratively establishes the best threshold for saccade detection for this
% trial. This is done either based on the eye velocity trace or based on
% the velocity after cross correlated with a saccade template.
% Given a threshold T, we take all data of the trial where the eye velocity
% as less than T and calculate its mean and std. We then establish a new
% threshold at V = mean + 6*std. We repeat this until the new threshold is
% less than 1°/s lower than the old threshold. The procedure is the same
% except that the procedure finished when the change in correlation
% threshold is less than 1% (.01).
% The thus established velocity threshold is close to optimal in that it
% allows us to detect as many saccades as possible given the noise in the
% data but does not detect too many due to possibly high noise in the data.
%
% !! Call this function with two parameters, unless you know what you are
% doing.

% prepare algorithm parameters
minFixSamples       = ceil(ETparams.fixation.minDur/1000 * ETparams.samplingFreq);
centralFixSamples   = ceil(ETparams.saccade.minDur /3000 * ETparams.samplingFreq);

% running (possibly partially) on xcorr output
if ETparams.data.qApplySaccadeTemplate
    [data.saccade.peakXCorrThreshold, meanData, stdData] = ...
        doOptimize(data.deg.velXCorr,ETparams.saccade.peakXCorrThreshold,.01, minFixSamples, centralFixSamples, nargin==2||qusecentralsample);
    
    if ETparams.saccade.qSaccadeTemplateRefine
        data.saccade.onsetXCorrThreshold = meanData + 3*stdData;
    end
end

% running on velocity trace, also if peaks are detected from the cross
% correlation trace.
% ETparams.saccade.qSaccadeTemplateRefine determines whether saccade
% onset/offset refinement is based on the velocity trace (false). If so,
% these are needed.
if ~ETparams.data.qApplySaccadeTemplate || ~ETparams.saccade.qSaccadeTemplateRefine
    [data.saccade.peakVelocityThreshold, meanData, stdData] = ...
        doOptimize(data.deg.vel,ETparams.saccade.peakVelocityThreshold,1, minFixSamples, centralFixSamples, nargin==2||qusecentralsample);
    
    data.saccade.onsetVelocityThreshold = meanData + 3*stdData;
end




function [peakThreshold, meanData, stdData] = doOptimize(data,initialThreshold,exitCriterion,minFixSamples, centralFixSamples, qUseCentralSamples)
% assign initial thresholds
peakThreshold = initialThreshold;
previousPeakDetectionThreshold = inf;

% iterate while we're gaining more than a 1° decrease in saccade peak
% velocity threshold
while previousPeakDetectionThreshold - peakThreshold > exitCriterion
    
    previousPeakDetectionThreshold = peakThreshold;
    
    % Find parts where the velocity is below the threshold, possible
    % fixation time (though this is still crude)
    qBelowThresh = data < peakThreshold;
    
    if qUseCentralSamples
        % We need to cut off the edges of the testing intervals and only
        % use parts of the data that are likely to belong to fixations or
        % the iteration will not converge to a lower threshold. So always
        % use this code path (just call this function with 2 arguments),
        % unless you want to see it for yourself. This is not just done to
        % speed up convergence.
        % NB: although this does not match how the algorithm is described
        % in Nyström & Holmqvist, 2010, it does match the code they made
        % available. As explained above, the more simple method they
        % described in their paper does not converge to lower thresholds
        % (in the minimal testing I did at least), while this code-path
        % appears to be robust.
        
        % get bounds of these detected peaks
        [threshon,threshoff] = bool2bounds(qBelowThresh);
        
        % throw out intervals that are too short and therefore unlikely to
        % be fixations
        qLongEnough = threshoff-threshon >= minFixSamples;
        threshon    = threshon (qLongEnough);
        threshoff   = threshoff(qLongEnough);
        
        % shrink them as done in Nystrom's version, to make sure we don't
        % catch the data that is still during the saccade
        threshon    = threshon  + centralFixSamples;
        threshoff   = threshoff - centralFixSamples;
        
        % convert to data selection indices
        idx         = bounds2ind(threshon,threshoff);
        
        % get mean and std of this data
        meanData    = nanmean(data(idx));
        stdData     = nanstd (data(idx));
    else
        meanData    = nanmean(data(qBelowThresh));
        stdData     = nanstd (data(qBelowThresh));
    end
    
    % calculate new threshold
    peakThreshold   = meanData + 6*stdData;
end