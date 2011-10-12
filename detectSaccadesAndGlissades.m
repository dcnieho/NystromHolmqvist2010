function data = detectSaccadesAndGlissades(data,ETparams)
% Detects start and end by criteria on either the eye velocity or the xcorr
% response between the velocity trace and a saccade template

% select parameter and data to work with
if ETparams.data.qApplySaccadeTemplate && ETparams.saccade.qSaccadeTemplateRefine
    % run full detection algorithm from the cross correlation trace
    field_peak  = 'peakXCorrThreshold';
    field_onset = 'onsetXCorrThreshold';
    field_offset= 'offsetXCorrThreshold';
    vel         = data.deg.velXCorr;
else
    % if ETparams.data.qApplySaccadeTemplate==true, then below just peaks
    % are detected based on xcorr responses, refinement is done from the
    % velocity trace (recommended). The detrended velocity trace is used
    % iff it is available and if the saccade template is not used
    field_peak  = 'peakVelocityThreshold';
    field_onset = 'onsetVelocityThreshold';
    field_offset= 'offsetVelocityThreshold';
    vel         = data.deg.vel;
end


%%% find where velocity data is above threshold
if ETparams.data.qApplySaccadeTemplate
    % find peaks from xcorr responses
    qAboveThresh    = data.deg.velXCorr > data.saccade.peakXCorrThreshold;
else
    qAboveThresh    = data.deg.vel      > data.saccade.peakVelocityThreshold;
end
[sacon,sacoff]  = bool2bounds(qAboveThresh);

%%% make vector that will contain true where saccades are detected
% this is used as some part so the algorithm (e.g. offset velocity
% threshold) needs to analyze data that is not during saccades or
% glissades.
qSacGlisOrNan   = isnan(vel);


% If no saccades are detected, return
if isempty(sacon)
    fprintf('no saccades\n');
    return;
end

% % NASA style saccade detection
% sacon   = sacon-10;
% sacoff  = sacoff+10;
% sacon(sacon<1)              = 1;
% sacoff(sacoff>length(vel))  = length(vel);

%%% prepare algorithm parameters
minSacSamples           = ceil(ETparams.saccade.minDur/1000                 * ETparams.samplingFreq);
maxGlissadeSamples      = ceil(ETparams.glissade.maxDur/1000                * ETparams.samplingFreq);
glissadeWindowSamples   = ceil(ETparams.glissade.searchWindow/1000          * ETparams.samplingFreq);
localNoiseWindowSamples = ceil(ETparams.saccade.localNoiseWindowLength/1000 * ETparams.samplingFreq);
% If the peak consists of =< minPeakSamples consequtive samples, it it
% probably noise (1/6 of the min saccade duration)
minPeakSamples          = ceil(ETparams.saccade.minDur/6000                 * ETparams.samplingFreq);
% if saccade is followed by another saccade by less than SacMergeWindow (in
% ms), they'll be merged
SacMergeWindowSamp      = ceil(ETparams.saccade.mergeWindow./1000           * ETparams.samplingFreq);


%%% Process one velocity peak at the time.
% Keep a counter here of how many peaks from sacon we have processed
% already. We need to use a while loop instead of a for-loop as the length
% of the saccade vector will change as we delete from it
kk = 1;

% for storing saccade end thresholds, which are adapted for each saccade
% based on noise before the saccade (not actually needed here, just put it
% here for documentation purposes).
saccadeOffsetTreshold  = [];

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
    while i > 1 && ...                                          % make sure we don't run out of the data
          (isnan(vel(i)) || ...                                 % and that we ignore nan data
           vel(i) > data.saccade.(field_onset) || ...           % keep searching until below saccadeOnsetVelocityTreshold
           diff(vel(i-[0:1])) < 0)                              % and acceleration is negative (we need to take this derivative locally as our acceleration signal is absolute)
        i = i-1;
    end
    sacon(kk) = i;                                              % velocity minimum is last sample before "acceleration" sign change
    
    % Calculate local noise during 'fixation' before the saccade start (the
    % adaptive part), excluding data that is during known saccade time or
    % missing.
    % This assumes the saccade is preceded the same eye velocity and noise
    % level as after the saccade - its the best we can do.
    % starting from the already refined saccade beginning we get
    % minFixSamples samples before the saccade, or as many as we can get,
    % that are not nan and not during saccade.
    idx = find(~qSacGlisOrNan(1:sacon(kk)-1),localNoiseWindowSamples,'last');
    localVelNoise = vel(idx);
    localVelNoise = mean(localVelNoise) + 3*std(localVelNoise);
        
    % Check whether the local velocity noise exceeds the peak velocity
    % threshold
    if ~isnan(localVelNoise) && localVelNoise < data.saccade.(field_peak)
        saccadeOffsetTreshold(kk) = localVelNoise*0.3 + data.saccade.(field_onset)*0.7; % 30% local + 70% global
    else
        saccadeOffsetTreshold(kk) = data.saccade.(field_onset);
    end
    
    % Detect saccade end. Walk forward from detected saccade start to find
    % where the velocity is below the saccadeOffsetTreshold and
    % where the acceleration is positive (which indicates a local minimum
    % in the velocity function)
    i = sacoff(kk);
    while i < length(vel) && ...                                % make sure we don't run out of the data
          (isnan(vel(i)) || ...                                 % and that we ignore nan data
           vel(i) > saccadeOffsetTreshold(kk) || ...            % keep searching until below saccadeOffsetTreshold
           diff(vel(i+[0:1])) < 0)                              % and acceleration is positive (we need to take this derivative locally as our acceleration signal is absolute)
        i = i+1;
    end
    sacoff(kk) = i;
    
    % now, delete all the next saccades that are enclosed by this potential
    % saccade, i.e., delete any saccade whose unrefined end of is before
    % the end of this saccade as they would converge to the same saccade
    % interval. We do this now befor the final checks below as if the
    % current saccade is deleted by this check, any later saccade that will
    % converge to the same interval would be deleted as well.
    while kk+1<=length(sacoff) &&...                            % make sure we don't run out of the data
          sacoff(kk+1) <= sacoff(kk)
        sacon (kk+1) = [];
        sacoff(kk+1) = [];
        continue;
    end
              
    % If the saccade contains NaN samples, delete it
    if ~ETparams.saccade.allowNaN && any(isnan(vel(sacon(kk):sacoff(kk))))
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
    % flag it in the trace
    qSacGlisOrNan(sacon(kk):sacoff(kk)) = true;
    
    
    %----------------------------------------------------------------------  
    % DETECT GLISSADE
    %----------------------------------------------------------------------
    
    if ETparams.glissade.qDetect
        % store glissades found in this interval, marked for post processing
        foundGlissadeOff            = [];

        % Search only for glissade peaks in a window after the saccade end
        % (both low and high velocity)
        % this window is also used below to assign the glissade type
        glissadeWindow = sacoff(kk) : min(sacoff(kk) + glissadeWindowSamples, length(vel)-1);

        % Detect glissade (low velocity criteria -- the adapted saccade offset
        % criterion)
        qGlissadePeak = vel(glissadeWindow) >= saccadeOffsetTreshold(kk);

        % Detect only 'complete' peaks (those with a beginning and an end)
        [~,potend] = bool2bounds(qGlissadePeak);
        % delete end of data that if it is detected -> glissade doesn't end
        % before end of window
        potend(potend==length(glissadeWindow)) = [];

        if ~isempty(potend)
            % found potential glissade, store it
            foundGlissadeOff  = sacoff(kk)+potend(end);                     % glissade start is saccade end
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
                        calcAmplitudeFick(data.deg.Azi(sacoff(kk)),data.deg.Ele(sacoff(kk)) , data.deg.Azi(sacoff(kk+idx)),data.deg.Ele(sacoff(kk+idx)))...
                      ) ...
                   < ...
                   max(...
                        calcAmplitudeFick(data.deg.Azi( sacon(kk)),data.deg.Ele( sacon(kk)) , data.deg.Azi(sacoff(kk)    ),data.deg.Ele(sacoff(kk)    ))...
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
            while foundGlissadeOff < length(vel) && ...                             % make sure we don't run out of the data
                  (isnan(vel(i)) || ...                                             % and that we ignore nan data
                   vel(foundGlissadeOff) > saccadeOffsetTreshold(kk) || ...         % keep searching until below saccadeOffsetTreshold
                   diff(vel(foundGlissadeOff+[0:1])) < 0)                           % and acceleration is positive (we need to take this derivative locally as our acceleration signal is absolute)
                foundGlissadeOff = foundGlissadeOff+1;
            end

            % See if glissade is valid according to our criteria
            % do not allow glissade duration > 80 ms AND
            % the glissade should not contain any NaN samples
            if foundGlissadeOff-sacoff(kk)+1 <= maxGlissadeSamples &&...
               (~ETparams.glissade.allowNaN && ~any(isnan(vel(sacoff(kk):foundGlissadeOff))))
                % store the glissade
                glissadeon  = [glissadeon   sacoff(kk)];
                glissadeoff = [glissadeoff  foundGlissadeOff];
                % determine type (use window of glissadeWindowSamples or less
                % if glissade is shorter) and store
                glissadeVelWindow = sacoff(kk) : min(sacoff(kk) + glissadeWindowSamples, foundGlissadeOff); % this cant go outside of data as foundGlissadeOff is always valid at this point
                if max(vel(glissadeVelWindow)) > data.saccade.(field_peak)
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
                
                % flag it in the trace
                qSacGlisOrNan(glissadeon(end):glissadeoff(end)) = true;
            end
        end
    end
    
    % increase counter, process next peak
    kk = kk+1;
end

% now deal with saccades that are too close together, fuse two saccades
% with little time between them
SacMergeWindowSamp

%%% output
data.saccade .on                = sacon;
data.saccade .off               = sacoff;
data.saccade .(field_offset)    = saccadeOffsetTreshold;
if ETparams.glissade.qDetect
    data.glissade.on            = glissadeon;
    data.glissade.off           = glissadeoff;
    data.glissade.type          = glissadetype;
end