function data = detectSaccadesAndGlissades(data,ETparams)
% Detects start and end by velocity criteria

%%% find where velocity data is above threshold
qAboveThresh    = data.deg.vel > data.saccade.peakVelocityThreshold;
[sacon,sacoff]  = bool2bounds(qAboveThresh);


% If no saccades are detected, return
if isempty(sacon)
    return;
end

%%% prepare algorithm parameters
% If the peak consists of =< minPeakSamples consequtive samples, it it
% probably noise (1/6 of the min saccade duration)
minPeakSamples          = ceil(ETparams.saccade.minDur/6000        * ETparams.samplingFreq);
minSacSamples           = ceil(ETparams.saccade.minDur/1000        * ETparams.samplingFreq);
minFixSamples           = ceil(ETparams.fixation.minDur/1000       * ETparams.samplingFreq);
maxGlissadeSamples      = ceil(ETparams.glissade.maxDur/1000       * ETparams.samplingFreq);
glissadeWindowSamples   = ceil(ETparams.glissade.searchWindow/1000 * ETparams.samplingFreq);

%%% Process one velocity peak at the time.
% Keep a counter here of how many peaks from sacon we have processed
% already. We need to use a while loop instead of a for-loop as the length
% of the saccade vector will change as we delete from it
kk = 1;

% for storing saccade end velocity threshold, which is computed adaptively
% for each saccade based on noise before the saccade (not actually needed
% here, just put it here for documentation purposes).
saccadeOffsetVelocityTreshold  = [];

% for storing glissade indices
glissadeon  = [];
glissadeoff = [];
glissadetype= [];   % 1 is low velocity, 2 is high velocity

while kk <= length(sacon)
    
    %----------------------------------------------------------------------  
    % Check the saccade peak samples to eliminate noise
    %----------------------------------------------------------------------
    
    % If the peak consists of =< minPeakSamples consequtive samples, it it
    % probably noise (1/6 of the min saccade duration), delete it
    if sacoff(kk)-sacon(kk) <= minPeakSamples
        sacon (kk) = [];
        sacoff(kk) = [];
        continue;
    end
    
    %----------------------------------------------------------------------
    % DETECT SACCADE - refine saccade beginnings and ends
    %----------------------------------------------------------------------       
    
    % Detect saccade start. Walk back from detected saccade start to find
    % where the velocity is below the saccadeVelocityThreshold (mean+3*std)
    % and where the acceleration is negative (which indicates a local
    % minimum in the velocity function).
    i = sacon(kk);
    while i > 1 && ...                                                  % make sure we don't run out of the data
          (isnan(data.deg.vel(i)) || ...                                % and that we ignore nan data
           data.deg.vel(i) > data.saccade.onsetVelocityTreshold || ...  % keep searching until below saccadeOnsetVelocityTreshold
           diff(data.deg.vel(i-[0:1])) < 0)                             % and acceleration is negative (we need to take this derivative locally as our acceleration signal is absolute)
        i = i-1;
    end
    sacon(kk) = i;%+1;                                                    % velocity minimum is last sample before "acceleration" sign change
    
    % Calculate local noise during 'fixation' before the saccade start (the
    % adaptive part)
    % this assumes the saccade is preceded by fixation (or at least by the
    % same eye velocity and noise level as after the saccade - its the best
    % we can do).
    localVelNoise = data.deg.vel(max(1,sacon(kk) - minFixSamples) : sacon(kk));
    localVelNoise = mean(localVelNoise) + 3*std(localVelNoise);
        
    % Check whether the local vel. noise exceeds the peak vel. threshold.
    if ~isnan(localVelNoise) && localVelNoise < data.saccade.peakVelocityThreshold
        saccadeOffsetVelocityTreshold(kk) = localVelNoise*0.3 + data.saccade.onsetVelocityTreshold*0.7; % 30% local + 70% global
    else
        saccadeOffsetVelocityTreshold(kk) = data.saccade.onsetVelocityTreshold;
    end
    
    % Detect saccade end. Walk forward from detected saccade start to find
    % where the velocity is below the saccadeOffsetVelocityTreshold and
    % where the acceleration is positive (which indicates a local minimum
    % in the velocity function)
    i = sacoff(kk);
    while i < length(data.deg.vel) && ...                               % make sure we don't run out of the data
          (isnan(data.deg.vel(i)) || ...                                % and that we ignore nan data
           data.deg.vel(i) > saccadeOffsetVelocityTreshold(kk) || ...   % keep searching until below saccadeOffsetVelocityTreshold
           diff(data.deg.vel(i+[0:1])) < 0)                             % and acceleration is positive (we need to take this derivative locally as our acceleration signal is absolute)
        i = i+1;
    end
    sacoff(kk) = i;
    
    % now, delete all the next saccades that are enclosed by this potential
    % saccade, i.e., delete any saccade whose unrefined end of is before
    % the end of this saccade as they would converge to the same saccade
    % interval. We do this now befor the final checks below as if the
    % current saccade is deleted by this check, any later saccade that will
    % converge to the same interval would be deleted as well.
    while kk+1<=length(sacoff) &&...                                    % make sure we don't run out of the data
          sacoff(kk+1) <= sacoff(kk)
        sacon (kk+1) = [];
        sacoff(kk+1) = [];
        continue;
    end
              
    % If the saccade contains NaN samples, delete it
    if ~ETparams.saccade.allowNaN && any(isnan(data.deg.vel(sacon(kk):sacoff(kk))))
        sacon (kk) = [];
        sacoff(kk) = [];
        continue;
    end
        
    % Make sure the saccade duration exceeds the minimum duration or delete
    % it
    if sacoff(kk)-sacon(kk)+1 < minSacSamples
        sacon (kk) = [];
        sacoff(kk) = [];
        continue;
    end
    
    %%%%
    % Done. All the above criteria are fulfilled, we've got a saccade.
    %%%%
    
    
    
    %----------------------------------------------------------------------  
    % DETECT GLISSADE
    %----------------------------------------------------------------------
    
    if ETparams.glissade.qDetect
        % store glissades found in this interval, marked for post processing
        foundGlissadeOff            = [];

        % Search only for glissade peaks in a window after the saccade end
        % (both low and high velocity)
        % this window is also used below to assign the glissade type
        glissadeWindow = sacoff(kk) : min(sacoff(kk) + glissadeWindowSamples, length(data.deg.vel)-1);

        % Detect glissade (low velocity criteria -- the adapted saccade offset
        % criterion)
        qGlissadePeak = data.deg.vel(glissadeWindow) >= saccadeOffsetVelocityTreshold(kk);

        % Detect only 'complete' peaks (those with a beginning and an end)
        [~,potend] = bool2bounds(qGlissadePeak);
        % delete end of data that if it is detected -> glissade doesn't end
        % before end of window
        potend(potend==length(glissadeWindow)) = [];

        if ~isempty(potend)
            % found potential glissade, store it
            foundGlissadeOff  = sacoff(kk)+potend(end);     % glissade start is saccade end
        end    

        % Detect glissade (high velocity criteria). These have already been
        % picked up by saccade detection, see if any have been
        % A high velocity glissade is detected only if glissade velocity goes
        % below the saccadePeakVelocityThreshold within the search window
        if kk+1<=length(sacoff)                                             % make sure we don't run out of data
            qHighGlisOff = sacoff(kk+1:end) <= glissadeWindow(end);         % see if any saccade offsets happens within our glissade detection window
            if any(qHighGlisOff)
                idx = find(qHighGlisOff,1,'last');

                % additional checks to see if this could be a potential
                % glissade:
                % - The amplitude of this potential glissade has to be lower
                %   than that of the saccade we're currently processing
                % - Also, if we've already detected a low velocity glissade,
                %   only replace with this high velocity glissade if detected
                %   end is later
                if max(...  % amplitude
                        calcAmplitudeFick(data.deg.X(sacoff(kk)),data.deg.Y(sacoff(kk)) , data.deg.X(sacoff(kk+idx)),data.deg.Y(sacoff(kk+idx)))...
                      ) ...
                   < ...
                   max(...
                        calcAmplitudeFick(data.deg.X( sacon(kk)),data.deg.Y( sacon(kk)) , data.deg.X(sacoff(kk)    ),data.deg.Y(sacoff(kk)    ))...
                      ) && ...
                   ((~isempty(foundGlissadeOff) && sacoff(kk+idx)>foundGlissadeOff) || isempty(foundGlissadeOff))   % end later then low velocity glissade
                    % found potential glissade, store it
                    foundGlissadeOff  = sacoff(kk+idx);
                end
            end
        end

        % If potential glissade detected
        if ~isempty(foundGlissadeOff)
            % Detect end. Walk forward from detected saccade start to find
            % where the acceleration is positive (which indicates a local
            % minimum in the velocity function)
            while foundGlissadeOff < length(data.deg.vel) && ...                                % make sure we don't run out of the data
                  (isnan(data.deg.vel(i)) || ...                                                % and that we ignore nan data
                   data.deg.vel(foundGlissadeOff) > saccadeOffsetVelocityTreshold(kk) || ...    % keep searching until below saccadeOffsetVelocityTreshold
                   diff(data.deg.vel(foundGlissadeOff+[0:1])) < 0)                              % and acceleration is positive (we need to take this derivative locally as our acceleration signal is absolute)
                foundGlissadeOff = foundGlissadeOff+1;
            end

            % See if glissade is valid according to our criteria
            % do not allow glissade duration > 80 ms AND
            % the glissade should not contain any NaN samples
            if foundGlissadeOff-sacoff(kk)+1 <= maxGlissadeSamples &&...
               (~ETparams.glissade.allowNaN && ~any(isnan(data.deg.vel(sacoff(kk):foundGlissadeOff))))
                % store the glissade
                glissadeon  = [glissadeon   sacoff(kk)];
                glissadeoff = [glissadeoff  foundGlissadeOff];
                % determine type (use window of glissadeWindowSamples or less
                % if glissade is shorter) and store
                glissadeVelWindow = sacoff(kk) : min(sacoff(kk) + glissadeWindowSamples, foundGlissadeOff); % this cant go outside of data as foundGlissadeOff is always valid at this point
                if max(data.deg.vel(glissadeVelWindow)) > data.saccade.peakVelocityThreshold
                    glissadetype = [glissadetype 2];    % high velocity glissade detected ('Strong glissade detection criteria')
                else
                    glissadetype = [glissadetype 1];    % low  velocity glissade detected ('Weak   glissade detection criteria')
                end
                % now, remove all next items in saccade vector that end before
                % end of this glissade. This also deletes the next saccade if
                % above we detected it to be a high velocity glissade
                while kk+1<=length(sacoff) &&...                            % make sure we don't run out of the data
                      sacoff(kk+1) <= glissadeoff(end)
                    sacon (kk+1) = [];
                    sacoff(kk+1) = [];
                    continue;
                end
            end
        end
    end
    
    % increase counter, process next peak
    kk = kk+1;
end

%%% output
data.saccade .on                        = sacon;
data.saccade .off                       = sacoff;
data.saccade .offsetVelocityTreshold    = saccadeOffsetVelocityTreshold;
data.glissade.on                        = glissadeon;
data.glissade.off                       = glissadeoff;
data.glissade.type                      = glissadetype;