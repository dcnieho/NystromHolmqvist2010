function detectSaccades(i,j)
% Detects start and end by velocity criteria

global ETparams

V = ETparams.data(i,j).vel;
A = ETparams.data(i,j).acc;
len = length(V);

% Preallocate memory
velLabeled = bwlabel(ETparams.data(i,j).velPeakIdx);
ETparams.saccadeIdx(i,j).Idx = zeros(1,len);    % Saccade index
ETparams.glissadeIdx(i,j).Idx = zeros(1,len);  % Glissade index

% If no saccades are detected, return
if isempty(velLabeled);
    return
end

% Process one velocity peak at the time
kk = 1;

for k = 1:max(velLabeled)
    
    %----------------------------------------------------------------------  
    % Check the saccade peak samples
    %----------------------------------------------------------------------       
    % The samples related to the current saccade
    peakIdx = find(velLabeled == k);
    
    % If the peak consists of =< minPeakSamples consequtive samples, it it probably
    % noise (1/6 or the min saccade duration)
    minPeakSamples = ceil(ETparams.minSaccadeDur/6*ETparams.samplingFreq); 
    if length(peakIdx) <= minPeakSamples, continue, end
    
    % Check whether this peak is already included in the previous saccade
    % (can be like this for glissades)
    if kk > 1
        if ~isempty(intersect(peakIdx,[find(ETparams.saccadeIdx(i,j).Idx) find(ETparams.glissadeIdx(i,j).Idx)]))
            continue
        end       
    end
       
    %----------------------------------------------------------------------
    % DETECT SACCADE
    %----------------------------------------------------------------------       
    
    % Detect saccade start.  AND acc <= 0
    saccadeStartIdx = find(V(peakIdx(1):-1:1) <= ETparams.data(i,j).saccadeVelocityTreshold &...% vel <= global vel threshold
                                                 [diff(V(peakIdx(1):-1:1)) 0] >= 0);          % acc <= 0
    if isempty(saccadeStartIdx), continue, end
    saccadeStartIdx = peakIdx(1) - saccadeStartIdx(1) + 1;
    
    % Calculate local fixation noise (the adaptive part)
    localVelNoise = V(saccadeStartIdx:-1: max(1,ceil(saccadeStartIdx - ETparams.minFixDur*ETparams.samplingFreq)));
    localVelNoise = mean(localVelNoise) + 3*std(localVelNoise);
    localsaccadeVelocityTreshold = localVelNoise*0.3 + ETparams.data(i,j).saccadeVelocityTreshold*0.7; % 30% local + 70% global
    
    % Check whether the local vel. noise exceeds the peak vel. threshold.
    if localVelNoise > ETparams.data(i,j).peakDetectionThreshold, continue, end
              
    % Detect end of saccade (without glissade)
    saccadeEndIdx = find(V(peakIdx(end):end) <= localsaccadeVelocityTreshold &...             % vel <= adaptive vel threshold
                                                  [diff(V(peakIdx(end):end)) 0] >= 0);        % acc <= 0
    
    if isempty(saccadeEndIdx), continue, end      
    saccadeEndIdx = peakIdx(end) + saccadeEndIdx(1) - 1;
    
    % If the saccade contains NaN samples, continue
    if any(ETparams.nanIdx(i,j).Idx(saccadeStartIdx:saccadeEndIdx)), continue, end
        
    % Make sure the saccade duration exceeds the minimum duration.
    saccadeLen = saccadeEndIdx - saccadeStartIdx;
    if saccadeLen/ETparams.samplingFreq < ETparams.minSaccadeDur
        continue    
    end
    
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

    %----------------------------------------------------------------------  
    % DETECT GLISSADE (ETparams.glissadeInfo(i,j,kk).type
    %----------------------------------------------------------------------   
    % Search only for glissade peaks in a window <= min fix duration after
    % the saccade end
    potentialGlissadeIdx = V(saccadeEndIdx: min(saccadeEndIdx + ETparams.minFixDur*ETparams.samplingFreq,len-1));
    
    % Detect glissade (low velocity criteria) 
    glissadePeakIdxW = potentialGlissadeIdx >= localsaccadeVelocityTreshold;    
    % Detect only 'complete' peaks (those with a beginning and an end)
    endIdx = find(abs(diff(glissadePeakIdxW)));   
    if length(endIdx)>1
        endIdx = endIdx(2:2:end); 
        glissadeEndWeakIdx = endIdx(end);
        nGlissadesWeak = length(endIdx);
    else
        glissadeEndWeakIdx = [];
    end
    
    % Detect glissade (high velocity criteria)
    glissadePeakIdxS = potentialGlissadeIdx >= ETparams.data(i,j).peakDetectionThreshold;
    glissadeEndStrongIdx = find(glissadePeakIdxS,1,'last');       
    nGlissadesStrong = length(unique(bwlabel(glissadePeakIdxS))) - 1;
    
    % Make sure that the saccade amplitude is larger than the glissade
    % amplitued, otherwise no glissade is detected.
    if max(potentialGlissadeIdx) > max(V(saccadeStartIdx:saccadeEndIdx)) 
        glissadeEndWeakIdx = [];
        glissadeEndStrongIdx = [];
    end
    
    % If no glissade detected
   if isempty(glissadeEndWeakIdx), 
        ETparams.glissadeInfo(i,j,kk).type = 0; 
        ETparams.glissadeInfo(i,j,kk).duration = 0;        
   % If glissade detected    
   else
       % Detect end.
        glissadeEndIdx = saccadeEndIdx + glissadeEndWeakIdx; 
        glissadeEndIdx = glissadeEndIdx + ...
                    find(diff(V(glissadeEndIdx:end)) >= 0,1,'first') - 1; 
        glissadeIdx = saccadeEndIdx:glissadeEndIdx;
        glissadeDuration = (length(glissadeIdx))/ETparams.samplingFreq;   
    
        
        % Do not allow glissade duration > 80 ms OR
        % If the glissade contains any NaN samples, continue
        if glissadeDuration > 2*ETparams.minFixDur |...
           any(ETparams.nanIdx(i,j).Idx(glissadeIdx))
            ETparams.glissadeInfo(i,j,kk).type = 0; 
            ETparams.glissadeInfo(i,j,kk).duration = 0;
            glissadeEndStrongIdx = [];
        else
            % Collect information about the glissade
            ETparams.glissadeInfo(i,j,kk).type = 1; % 'Weak glissade detection criteria'
            ETparams.glissadeIdx(i,j).Idx(glissadeIdx) = 1;
            ETparams.glissadeInfo(i,j,kk).duration = glissadeDuration;
            ETparams.glissadeInfo(i,j,kk).numberOf = [nGlissadesWeak, nGlissadesStrong];
        end
        
   end
    
   if ~isempty(glissadeEndStrongIdx)
        ETparams.glissadeInfo(i,j,kk).type = 2; % 'Strong glissade detection criteria'
   end                                                                                  

    kk = kk+1;
    
end

