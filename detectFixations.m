function data = detectFixations(data,ETparams)
%--------------------------------------------------------------------------
% Fixation detection
% Fixation are detected implicitly, as sections that are not saccade or
% glissade and also fit some other criteria
%--------------------------------------------------------------------------

%%% find data that is not saccade or glissade
[fixon,fixoff]  = findContiguousRegions(~(...
                    bounds2bool(data.saccade.on ,data.saccade.off ,length(data.deg.vel)) | ...
                    bounds2bool(data.glissade.on,data.glissade.off,length(data.deg.vel))   ...
                  ));

%%% prepare algorithm parameters
minFixSamples   = ceil(ETparams.fixation.minDur/1000 * ETparams.samplingFreq);

%%% Process the tentative fixations
% Keep a counter here of how many sections from fixon we have processed
% already. We need to use a while loop instead of a for-loop as the length
% of the fixation vector will change as we delete from it
kk = 1;

while kk <= length(fixon)
    % Check that the fixation duration exceeds the minimum duration
    % criteria. Delete if not
    if fixoff(kk)-fixon(kk) < minFixSamples
        fixon (kk) = [];
        fixoff(kk) = [];
        continue;
    end
    
    % Exclude section if any of the samples has a velocity > peak saccade
    % threshold, it cannot be a fixation (section somehow got deleted by
    % the saccade algorithm)
    if any(data.deg.vel(fixon(kk):fixoff(kk)) > data.saccade.peakVelocityThreshold)
        fixon (kk) = [];
        fixoff(kk) = [];
        continue;
    end
    
    % If the fixation contains NaN samples, delete it
    if any(isnan(data.deg.vel(fixon(kk):fixoff(kk))))
        fixon (kk) = [];
        fixoff(kk) = [];
        continue;
    end
    
    %%%%
    % Done. All the above criteria are fulfilled, we've got a fixation.
    %%%%
    
    % increase counter, process next section
    kk = kk+1;
end

%%% output
data.fixation.on    = fixon;
data.fixation.off   = fixoff;