function data = detectAndRemoveBlinks(data,ETparams)
% Detects blinks and other un-physiological eye movements. They are
% replaced with nan
% 
% Blinks are detect in one of two ways.
% 1. Velocity or acceleration trace shows speeds above what is
%    physiologically possible.
% 2. High-amplitude downward excursions of recorded gaze position usually
%    accompanied with outlier speeds. Furthermore, these have an
%    impossible downward peak-like shape in the position trace.


% prepare parameters
% if blink is followed by another blink by less than mergeWindow (in ms),
% they'll be merged
blinkMergeWindowSamples = ceil(ETparams.blink.mergeWindow./1000 * ETparams.samplingFreq);


% TODO: first, check NaNs are in the same places in all traces. Sanity data integrity check
fn = fieldnames(data.deg);
% all position fields should have same NaNs, ass should all velocity fields
% and all acceleration fields

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method 1: unphysiological eye movements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detect possible blinks, episodes where the eyes move too fast too be
% physiologically possible
qBlink =     data.deg.vel  > ETparams.blink.velocityThreshold |...
         abs(data.deg.acc) > ETparams.blink.accThreshold;

% find bounds of blinks or noise as detected above
[blinkon,blinkoff]      = bool2bounds(qBlink);

% Process one possible blink at a time, refine the bounds
% These unphysiological movements should already have been detected as
% saccades (they have very large velocity!). So use detected saccade on and
% offsets as blink on- and offsets.
sacon  = data.saccade.on;
sacoff = data.saccade.off;
qIsBlink = false(size(sacon));
for p = 1:length(blinkon)
    qEnclosed = sacon<=blinkon(p) & sacoff>=blinkoff(p);
    
    if any(qEnclosed)
        assert(sum(qEnclosed)==1)  % anything else would be ridiculous at this stage as no overlapping saccades should exist!
        blinkon(p)  = sacon (qEnclosed);
        blinkoff(p) = sacoff(qEnclosed);
        
        % mark saccade as blink
        qIsBlink(qEnclosed) = true;
        
        if 0
            % debug
            figure(200),clf
            plot(data.deg.Ele(blinkon(p):blinkoff(p)))
            title(sprintf('blink %d',p))
            pause
        end
    else
        % all blinks should have been detected as saccades in the previous
        % steps
        error('Not implemented/shouldn''t occur. Find blink start end in another way.')
    end
end

% multiple unphysiological segments might be enclosed by same saccade, in
% which case the same blink is detected multiple times. Remove duplicated
[~,idx] = unique(blinkon);
blinkon  = blinkon (idx);
blinkoff = blinkoff(idx);
assert(all(blinkoff>blinkon) && all(blinkoff(1:end-1)<blinkon(2:end)))  % sanity checks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method 2: downward peak-like structures in the vertical/elevation
%           position trace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% each saccade should be a roughly monotonic position shift. If instead teh
% eye ends up back where it started vertically (at least w.r.t. its maximum
% downward position during the saccade, we're dealing with another class of
% un-physiological movement that indicates a blink-like eye recording
% wobble (apparently during blink-like eye movements the Eyelink doesn't
% always emit missing samples!)
for p = 1:length(sacon)
    if qIsBlink(p)
        % already flagged as blink, skip
        continue;
    end
    
    % get vertical eye position during saccade
    eye_pos = data.deg.Ele(sacon(p):sacoff(p));
    
    % subtract position at begin so we have relative positions during
    % saccade
    eye_pos = eye_pos - eye_pos(1);

    % get max amplitude downward vertical position (downward is positive)
    eye_max = max(eye_pos);
    
    % get amplitude at end
    eye_end = mean(eye_pos(find(~isnan(eye_pos),10,'last')));
    
    % if:
    % 1. saccade is downward
    % 2. it is downward by at least 5 degrees
    % 3. final amplitude is less than 40% of max amplitude (the eye apparently turned around)
    % -> We're dealing with a blink
    if nanmean(eye_pos)>0 && ...    % 1
            eye_max > 5 && ...      % 2
            eye_end < eye_max*.4    % 3
        
        % add info about blink
        blinkon  = [blinkon  sacon(p) ];
        blinkoff = [blinkoff sacoff(p)];
        % mark saccade as blink
        qIsBlink(p) = true;
        
        if 0
            % debug
            figure(200),clf
            plot(eye_pos)
            title(sprintf('saccade %d',p))
            pause
        end
    end
end

% build information about blinks
[blinkon,idx]   = sort(blinkon);
data.blink.on   = blinkon;
data.blink.off  = blinkoff(idx);

% now remove saccades that were flagged as blink
% first remove their corresponding glissades, if any
if ETparams.glissade.qDetect
    qRemoveGlissade = ismember(data.glissade.on,data.saccade.off(qIsBlink));
    data.glissade   = removeElementFromStructFields(data.glissade,qRemoveGlissade);
end

% then deal with the saccades
data.saccade    = removeElementFromStructFields(data.saccade,qIsBlink);

% merge very close blinks.
data.blink      = mergeIntervals(data.blink, [], blinkMergeWindowSamples);

% replace with nan if wanted
if ETparams.blink.qReplaceWithNan
    % create boolean matrix given blink bounds
    qBlink = bounds2bool(data.blink.on+1,data.blink.off-1,length(data.deg.vel));    % remove one sample inwards as thats good for plotting and otherwise doesn't matter
    % remove data that is due to noise
    data.deg.vel(qBlink)    = nan;
    % TODO: other traces?
    
    % lastly, notify if more than 20% nan
    if sum(isnan(data.deg.vel))/length(data.deg.vel) > 0.20
        fprintf('Warning: This trial contains %.2f%% missing+blinks samples\n',sum(isnan(data.deg.vel))/length(data.deg.vel)*100);
        data.qNoiseTrial = true;
    else
        data.qNoiseTrial = false;
    end
end