function detectFixations(i,j)
%--------------------------------------------------------------------------
% Fixation detection
% Fixation are detected implicitly 
%--------------------------------------------------------------------------

global ETparams

possibleFixationIdx = ~(ETparams.saccadeIdx(i,j).Idx | ETparams.glissadeIdx(i,j).Idx);
fixLabeled = bwlabel(possibleFixationIdx);

% Process one fixation (or more precisely a period of gaze samples that might
% a fixation) at the time.
kk = 1;
ETparams.fixationIdx(i,j).Idx = zeros(1,length(ETparams.data(i,j).X));
for k = 1:max(fixLabeled)

    % The samples related to the current fixation
    fixIdx = find(fixLabeled == k);
    
    % Check that the fixation duration exceeds the minimum duration criteria. 
    if length(fixIdx)/ETparams.samplingFreq < ETparams.minFixDur
        continue    
    end
       
    % If any of the sample has a velocity > peak saccade threshold, it
    % cannot be a fixation (missed by the saccade algorithm)
    if any(ETparams.data(i,j).vel(fixIdx) > ETparams.data(i,j).peakDetectionThreshold)
        continue
    end
    
    % If the saccade contains NaN samples, continue
    if any(ETparams.nanIdx(i,j).Idx(fixIdx)), continue, end
    
    % If all the above criteria are fulfilled, label it as a fixation.
    ETparams.fixationIdx(i,j).Idx(fixIdx) = 1;
    
    % Calculate the position of the fixation
    ETparams.fixationInfo(i,j,kk).X = nanmean(ETparams.data(i,j).X(fixIdx));
    ETparams.fixationInfo(i,j,kk).Y = nanmean(ETparams.data(i,j).Y(fixIdx));

    % Collect information about the fixation
    fixationStartIdx = fixIdx(1);
    fixationEndIdx = fixIdx(end);

    ETparams.fixationInfo(i,j,kk).start = fixationStartIdx/ETparams.samplingFreq; % in ms
    ETparams.fixationInfo(i,j,kk).end = fixationEndIdx/ETparams.samplingFreq; % in ms
    ETparams.fixationInfo(i,j,kk).duration = ETparams.fixationInfo(i,j,kk).end - ETparams.fixationInfo(i,j,kk).start;
    
    kk = kk+1;
    
end

