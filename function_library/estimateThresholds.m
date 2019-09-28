function data = estimateThresholds(data,ETparams,useCentralSample)
% iteratively establishes the best threshold for saccade classification for
% this trial. This is done either based on the eye velocity trace or based
% on the output of cross correlation with a saccade template.
% Given a threshold T, we take all data of the trial where the eye velocity
% is less than T and calculate its mean and std. We then establish a new
% threshold at V = mean + 6*std. We repeat this until the new threshold is
% less than 1°/s lower than the old threshold. The procedure is the same
% when using xcorr output, except that the procedure finishes when the
% change in correlation threshold is less than .01.
% The thus established velocity threshold is close to optimal in that it
% allows us to classify as many saccades as possible given the noise in the
% data but does not generate many false alarms when there is high noise in
% the data.
%
% !! Call this function with two parameters, unless you know what you are
% doing.

% prepare algorithm parameters
minFixSamples       = ceil(50                      /1000 * ETparams.samplingFreq);  % min duration of possible fixation (non-saccade actually) for it to be included in this estimation process
centralFixSamples   = ceil(ETparams.saccade.minDur /3000 * ETparams.samplingFreq);

% saccade classification will run (possibly partially) on xcorr output,
% need thresholds from this trace
if ETparams.data.applySaccadeTemplate
    [data.saccade.peakXCorrThreshold, meanData, stdData] = ...
        doOptimize(data.deg.velXCorr,ETparams.saccade.peakXCorrThreshold,.01, ETparams.saccade.peakXCorrSD, minFixSamples, centralFixSamples, nargin==2||useCentralSample);
    
    if ETparams.saccade.useTemplateForRefine
        data.saccade.onsetXCorrThreshold = meanData + 3*stdData;
    end
end

% saccade classification will run on the velocity trace, possibly while
% peaks are classified from the cross correlation trace.
% ETparams.saccade.useTemplateForRefine determines whether saccade
% onset/offset refinement is based on the velocity trace (false). If so,
% these are needed.
if ~ETparams.data.applySaccadeTemplate || ~ETparams.saccade.useTemplateForRefine
    [data.saccade.peakVelocityThreshold, meanData, stdData] = ...
        doOptimize(data.deg.vel,ETparams.saccade.peakVelocityThreshold,1, ETparams.saccade.peakVelocitySD, minFixSamples, centralFixSamples, nargin==2||useCentralSample);
    
    data.saccade.onsetVelocityThreshold = meanData + 3*stdData;
end

% also do change of pupil size for blink classification
if bitand(ETparams.blink.classifyMode,uint8(1)) && isfield(data,'pupil') && isfield(data.pupil,'dsize')
    [data.blink.peakDSizeThreshold, meanData, stdData] = ...
        doOptimize(abs(data.pupil.dsize),ETparams.blink.dSizeThreshold,5, ETparams.blink.dSizeSD, minFixSamples, centralFixSamples, nargin==2||useCentralSample);

    data.blink.onsetDSizeThreshold = meanData + 3*stdData;
end



function [peakThreshold, meanData, stdData] = doOptimize(data,initialThreshold,exitCriterion, nStd, minFixSamples, centralFixSamples, qUseCentralSamples)
% assign initial thresholds
peakThreshold                       = initialThreshold;
previousPeakClassificationThreshold = inf;

% iterate while we're gaining more than a exitCriterion decrease in saccade peak
% velocity threshold
while previousPeakClassificationThreshold - peakThreshold > exitCriterion
    
    previousPeakClassificationThreshold = peakThreshold;
    
    % Find parts where the velocity is below the threshold, possible
    % data during fixation (though this is still crude)
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
        % available. As mentioned above, the more simple method they
        % described in their paper does not converge to lower thresholds
        % (in the minimal testing I did at least), while this code-path
        % appears to be robust.
        
        % get bounds of these found peaks
        [threshon,threshoff] = bool2bounds(qBelowThresh);
        
        % throw out intervals that are too short and therefore unlikely to
        % be fixations
        qLongEnough = threshoff-threshon >= minFixSamples;
        threshon    = threshon (qLongEnough);
        threshoff   = threshoff(qLongEnough);
        
        % shrink them as done in Nyström's version, to make sure we don't
        % catch the data that is still during the saccade
        threshon    = threshon  + centralFixSamples;
        threshoff   = threshoff - centralFixSamples;
        
        % convert to data selection indices
        idx         = bounds2ind(threshon,threshoff);
        
        % get mean and std of this data
        meanData    = mean(data(idx),'omitnan');
        stdData     = std (data(idx),'omitnan');
    else
        meanData    = mean(data(qBelowThresh),'omitnan');
        stdData     = std (data(qBelowThresh),'omitnan');
    end
    
    % calculate new threshold
    peakThreshold   = meanData + nStd*stdData;
    
    if isempty(idx)
        fprintf('couldn''t determine threshold\n');
    end
end