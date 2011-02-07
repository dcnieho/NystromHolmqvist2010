function data = detectFixations(data,ETparams)
%--------------------------------------------------------------------------
% Fixation detection
% Fixation are detected implicitly, as sections that are not saccade or
% glissade and also fit some other criteria
%--------------------------------------------------------------------------

%%% find data that is not saccade or glissade
[fixon,fixoff]  = findContiguousRegions(~(...
                    bounds2bool(data.saccade.on ,data.saccade.off ,length(data.vel)) | ...
                    bounds2bool(data.glissade.on,data.glissade.off,length(data.vel))   ...
                  ));

%%% prepare algorithm parameters
minFixSamples   = ceil(ETparams.minFixDurms/1000 * ETparams.samplingFreq);

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
    if any(data.vel(fixon(kk):fixoff(kk)) > data.peakDetectionThreshold)
        fixon (kk) = [];
        fixoff(kk) = [];
        continue;
    end
    
    % If the fixation contains NaN samples, delete it
    if any(isnan(data.vel(fixon(kk):fixoff(kk))))
        fixon (kk) = [];
        fixoff(kk) = [];
        continue;
    end
    
    %%%%
    % Done. All the above criteria are fulfilled, we've got a fixation.
    %%%%
    
    if 0
        % TODO: I want to put all this in some other step. Now we're just
        % establishing beginnings and ends in the data. This information
        % collection is logically a separate analysis step, and not always
        % needed
        
        % Calculate the position of the fixation
        ETparams.fixationInfo(i,j,kk).X = nanmean(ETparams.data(i,j).X(fixIdx));
        ETparams.fixationInfo(i,j,kk).Y = nanmean(ETparams.data(i,j).Y(fixIdx));
        
        % Collect information about the fixation
        fixationStartIdx = fixIdx(1);
        fixationEndIdx = fixIdx(end);
        
        ETparams.fixationInfo(i,j,kk).start = fixationStartIdx/ETparams.samplingFreq; % in ms
        ETparams.fixationInfo(i,j,kk).end = fixationEndIdx/ETparams.samplingFreq; % in ms
        ETparams.fixationInfo(i,j,kk).duration = ETparams.fixationInfo(i,j,kk).end - ETparams.fixationInfo(i,j,kk).start;
    end
    
    % increase counter, process next section
    kk = kk+1;
end

%%% output
data.fixation.on    = fixon;
data.fixation.off   = fixoff;