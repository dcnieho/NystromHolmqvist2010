function data = detectSaccades(data,ETparams)
% Detects start and end by velocity criteria

%%% find where velocity data is above threshold
qAboveThresh    = data.vel > data.peakDetectionThreshold;
[sacon,sacoff]  = findContiguousRegions(qAboveThresh);

% If no saccades are detected, return
if isempty(sacon)
    return;
end

%%% prepare algorithm parameters
% If the peak consists of =< minPeakSamples consequtive samples, it it
% probably noise (1/6 of the min saccade duration)
minPeakSamples          = ceil(ETparams.minSaccadeDurms/6000        * ETparams.samplingFreq);
minFixSamples           = ceil(ETparams.minFixDurms/1000            * ETparams.samplingFreq);
glissadeWindowSamples   = ceil(ETparams.glissadeSearchWindowms/1000 * ETparams.samplingFreq);

%%% Process one velocity peak at the time.
% Keep a counter here of how many peaks from sacon we have processed
% already as we might sometimes process two in the loop (in case of high
% velocity glissades)
kk = 0;

% keep boolean that marks which need to be deleted. Due to the possibility
% that a glissade follows a saccade, we can't run this backwards. Thus need
% to delete the crap afterwards
idxNeedDelete = [];

% for storing glissade indices
glissadeon  = [];
glissadeoff = [];
glissadetype= [];   % 1 is low velocity, 2 is high velocity

while kk+1 <= length(sacon)
    
    % increase counter, process next peak
    kk = kk+1;
    
    %----------------------------------------------------------------------  
    % Check the saccade peak samples to eliminate noise
    %----------------------------------------------------------------------
    
    % If the peak consists of =< minPeakSamples consequtive samples, it it
    % probably noise (1/6 of the min saccade duration), delete it
    if sacoff(kk)-sacon(kk) <= minPeakSamples
        idxNeedDelete = [idxNeedDelete kk];
        continue;
    end
    
    
    %----------------------------------------------------------------------
    % DETECT SACCADE
    %----------------------------------------------------------------------       
    
    % Detect saccade start. Walk back from detected saccade start to find
    % where the velocity is below the saccadeVelocityThreshold (mean+3*std)
    % and where the acceleration is negative (which indicates a local
    % minimum in the velocity function).
    i = sacon(kk);
    while i > 1 && ...                                          % make sure we don't run out of the data
          (data.vel(i) > data.saccadeVelocityTreshold || ...    % keep searching until below saccadeVelocityTreshold
           data.vel(i+1)-data.vel(i) > 0)                       % and acceleration is negative (we need to take this derivative locally as our acceleration signal is absolute)
        i = i-1;
    end
    sacon(kk) = i+1;    % velocity minimum is last sample before "acceleration" sign change
          
    % Calculate local fixation noise before the saccade start (the adaptive
    % part)
    % this assumes the saccade is preceded by fixation (or at least by the
    % same eye velocity and noise level as after the saccade - its the best
    % we can do).
    localVelNoise = data.vel(max(1,sacon(kk) - minFixSamples) : sacon(kk));
    localVelNoise = mean(localVelNoise) + 3*std(localVelNoise);
        
    % Check whether the local vel. noise exceeds the peak vel. threshold.
    if localVelNoise > data.peakDetectionThreshold, continue, end
    if ~isnan(localVelNoise)
        localsaccadeVelocityTreshold = localVelNoise*0.3 + data.saccadeVelocityTreshold*0.7; % 30% local + 70% global
    else
        localsaccadeVelocityTreshold = data.saccadeVelocityTreshold;
    end
    
    % Detect saccade end. Walk forward from detected saccade start to find
    % where the velocity is below the localsaccadeVelocityTreshold and
    % where the acceleration is positive (which indicates a local minimum
    % in the velocity function)
    i = sacoff(kk);
    while i < length(data.vel)-1 && ...                         % make sure we don't run out of the data
          (data.vel(i) > localsaccadeVelocityTreshold || ...    % keep searching until below localsaccadeVelocityTreshold
           data.vel(i+1)-data.vel(i) < 0)                       % and acceleration is positive (we need to take this derivative locally as our acceleration signal is absolute)
        i = i+1;
    end
    sacoff(kk) = i;
              
    % If the saccade contains NaN samples, delete it
    if any(isnan(data.vel(sacon(kk):sacoff(kk))))
        idxNeedDelete = [idxNeedDelete kk];
        continue;
    end
        
    % Make sure the saccade duration exceeds the minimum duration or delete
    % it
    if (sacoff(kk)-sacon(kk))/ETparams.samplingFreq < ETparams.minSaccadeDur
        idxNeedDelete = [idxNeedDelete kk];
        continue;
    end
    
    if 0
        % TODO: I want to put all this in some other step. Now we're just
        % astablishing beginnings and ends in the data. This information
        % collection is logically a separate analysis step, and not always
        % needed
        
        % If all the above criteria are fulfilled, label it as a saccade.
        ETparams.saccadeIdx(i,j).Idx(saccadeStartIdx:saccadeEndIdx) = 1;
        ETparams.data(i,j,kk).localSaccadeVelocityTreshold = localsaccadeVelocityTreshold;
        
        % Collect information about the saccade
        ETparams.saccadeInfo(i,j,kk).start = saccadeStartIdx/ETparams.samplingFreq; % in ms
        ETparams.saccadeInfo(i,j,kk).end = saccadeEndIdx/ETparams.samplingFreq; % in ms
        ETparams.saccadeInfo(i,j,kk).duration = ETparams.saccadeInfo(i,j,kk).end - ETparams.saccadeInfo(i,j,kk).start;
        ETparams.saccadeInfo(i,j,kk).amplitude = sqrt(((ETparams.data(i,j).X(saccadeEndIdx)-...
            (ETparams.data(i,j).X(saccadeStartIdx))))^2 + ...
            ((ETparams.data(i,j).Y(saccadeEndIdx)-...
            (ETparams.data(i,j).Y(saccadeStartIdx))))^2   );
        ETparams.saccadeInfo(i,j,kk).peakVelocity = max(V(saccadeStartIdx:saccadeEndIdx));
        ETparams.saccadeInfo(i,j,kk).peakAcceleration = max(A(saccadeStartIdx:saccadeEndIdx));
    end

    
    %----------------------------------------------------------------------  
    % DETECT GLISSADE (ETparams.glissadeInfo(i,j,kk).type
    %----------------------------------------------------------------------
    
    % store glissades found in this interval, marked for post processing
    foundGlissadeOff  = [];
    foundGlissadeType = [];
    
    % Search only for glissade peaks in a window after the saccade end
    % (both low and high velocity)
    % this window is also used below to assign the glissade type
    glissadeWindow = sacoff(kk) : min(sacoff(kk) + glissadeWindowSamples, length(data.vel)-1);
    
    % Detect glissade (low velocity criteria -- the adapted saccade offset
    % criterion)
    qGlissadePeak = data.vel(glissadeWindow) >= localsaccadeVelocityTreshold;
    
    % Detect only 'complete' peaks (those with a beginning and an end)
    [~,potend]  = findContiguousRegions(qGlissadePeak);
    % delete end of data that if it is detected -> glissade doesn't end
    % before end of window
    potend(potend==length(glissadeWindow)) = [];
    
    if ~isempty(potend)
        % found potential glissade, store it
        foundGlissadeOff  = sacoff(kk)+potend(end);     % glissade start is saccade end
    end
    
    % Detect glissade (high velocity criteria). This should have already
    % been picked up by saccade detection, see if it has been
    % A high velocity glissade is detected only if glissade velocity goes
    % below the peakDetectionThreshold within the search window
    % Also, make sure that the saccade amplitude is larger than the
    % glissade amplitude, otherwise no glissade is detected. (This operates
    % on max velocity instead, but thats fine due to the mean sequence)
    if kk+1<=length(sacoff) &&...                                                   % make sure we don't run out of data
       sacoff(kk+1) > sacoff(kk) && sacoff(kk+1) <= glissadeWindow(end) && ...      % see if next saccade offset happens within our glissade detection window
       ((~isempty(foundGlissadeOff) && sacoff(kk+1)>foundGlissadeOff) || isempty(foundGlissadeOff)) && ...  % if we've already detected a low velocity glissade, only replace with this high velocity glissade if detected end is later
       max(data.vel(sacoff(kk):sacoff(kk+1))) < max(data.vel(sacon(kk):sacoff(kk))) % and see if the peak velocity during this second saccade is lower than during the one processed above
        % found potential glissade, store it
        foundGlissadeOff  = sacoff(kk+1);       % glissade start is saccade end
    end
    
    % If glissade detected
    if ~isempty(foundGlissadeOff)
        % Detect end. Walk forward from detected saccade start to find
        % where the acceleration is positive (which indicates a local
        % minimum in the velocity function)
        while foundGlissadeOff < length(data.vel)-1 && ...                      % make sure we don't run out of the data
             (data.vel(foundGlissadeOff) > localsaccadeVelocityTreshold || ...  % keep searching until below localsaccadeVelocityTreshold
              data.vel(foundGlissadeOff+1)-data.vel(foundGlissadeOff) < 0)      % and acceleration is positive (we need to take this derivative locally as our acceleration signal is absolute)
            foundGlissadeOff = foundGlissadeOff+1;
        end
        
        glissadeDuration = (foundGlissadeOff-sacoff(kk))/ETparams.samplingFreq;
        
        % Do not allow glissade duration > 80 ms AND
        % the glissade should not contain any NaN samples
        if glissadeDuration <= 2*ETparams.minFixDur ||...
           ~any(isnan(data.vel(sacoff(kk):foundGlissadeOff)))
            % store the glissade
            glissadeon  = [glissadeon   sacoff(kk)];
            glissadeoff = [glissadeoff  foundGlissadeOff];
            % determine type (use window of glissadeWindowSamples or less
            % if glissade is shorter) and store
            glissadeVelWindow = sacoff(kk) : min(sacoff(kk) + glissadeWindowSamples, foundGlissadeOff); % this cant go outside of data as foundGlissadeOff is always valid at this point
            if max(data.vel(glissadeVelWindow)) > data.peakDetectionThreshold
                glissadetype = [glissadetype 2];    % high velocity glissade detected ('Strong glissade detection criteria')
            else
                glissadetype = [glissadetype 1];    % low  velocity glissade detected ('Weak   glissade detection criteria')
            end
            if glissadetype(end)==2
                % this glissade is also indicated in the saccade vector,
                % the next one from the current saccade to be sure
                % increase kk by one to skip it on the next loop
                % also, mark it for deletion from the saccades
                kk = kk+1;
                idxNeedDelete = [idxNeedDelete kk]; % we just increased kk, this is correct
            end
        end
    end
end

% post-process: delete saccades marked for deletion
sacon (idxNeedDelete) = [];
sacoff(idxNeedDelete) = [];

% output
data.sacon          = sacon;
data.sacoff         = sacoff;
data.glissadeon     = glissadeon;
data.glissadeoff    = glissadeoff;
data.glissadetype   = glissadetype;