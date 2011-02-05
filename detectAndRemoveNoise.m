function detectAndRemoveNoise(i,j)
% Detects and removes un-physiological movement (which derives from noise
% and blinks)
global ETparams

ETparams.nanIdx(i,j).Idx = zeros(1,length(ETparams.data(i,j).X));

V = ETparams.data(i,j).vel;
V_threshold = median(ETparams.data(i,j).vel)*2;

% Detect possible blinks and noise (where XY-coords are 0  or if the eyes move too fast)
blinkIdx = (ETparams.data(i,j).X <= 0 & ETparams.data(i,j).Y <= 0) |...
        ETparams.data(i,j).vel > ETparams.blinkVelocityThreshold |...
        abs(ETparams.data(i,j).acc) > ETparams.blinkAccThreshold;

% Set possible blink and noise index to '1'
ETparams.nanIdx(i,j).Idx(blinkIdx) = 1;   

% Label blinks or noise
blinkLabeled = bwlabel(blinkIdx);

% Process one blink or noise period at the time
for k = 1:max(blinkLabeled)

    % The samples related to the current event
    b = find(blinkLabeled == k);
      
    % Go back in time to see where the blink (noise) started
    sEventIdx = find(V(b(1):-1:1) <= V_threshold);
    if isempty(sEventIdx), continue, end
    sEventIdx = b(1) - sEventIdx(1) + 1;
    ETparams.nanIdx(i,j).Idx(sEventIdx:b(1)) = 1;      
    
    % Go forward in time to see where the blink (noise) started    
    eEventIdx = find(V(b(end):end) <= V_threshold);
    if isempty(eEventIdx), continue, end    
    eEventIdx = (b(end) + eEventIdx(1) - 1);
    ETparams.nanIdx(i,j).Idx(b(end):eEventIdx) = 1;
    
end

temp_idx = find(ETparams.nanIdx(i,j).Idx);
if temp_idx/length(V) > 0.20
    disp('Warning: This trial contains > 20 % noise+blinks samples')
    ETparams.data(i,j).NoiseTrial = 0;
else
    ETparams.data(i,j).NoiseTrial = 1;
end
ETparams.data(i,j).vel(temp_idx) = nan;
ETparams.data(i,j).acc(temp_idx) = nan;
ETparams.data(i,j).X(temp_idx) = nan;
ETparams.data(i,j).Y(temp_idx) = nan;


