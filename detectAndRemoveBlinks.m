function data = detectAndRemoveBlinks(data,ETparams)
% Detects blinks and other un-physiological eye movements. They are
% replaced with nan
% 
% Blinks are detect in one of two ways.
% 1. Velocity or acceleration trace shows speeds above what is
%    physiologically possible.
% 2. With the Eyelink, blinks are usually evident as high-amplitude
%    downward excursions of recorded gaze position, usually accompanied
%    with outlier speeds. Furthermore, these have an impossible downward
%    peak-like shape in the position trace.
%
% All of this code is kind of a kludge that happens to work for our
% dataset. Redo properly once....


% prepare parameters
% if blink is followed by another blink by less than mergeWindow (in ms),
% they'll be merged
blinkMergeWindowSamples = ceil(ETparams.blink.mergeWindow./1000 * ETparams.samplingFreq);


% First a sanity data integrity check: check NaN positionss are in the same
% places in all traces, if not I must have introduced some bug somewhere.
% Not that its important for the correct functioning of this function by
% the way...
velFieldsDeg = checkNansSameForAllFields(data,'deg');

% if we have smoothed eye position in pixels and its derivatives, check as
% well
if ETparams.data.qAlsoStoreandDiffPixels
    velFieldsPix = checkNansSameForAllFields(data,'pix');
    % and check correspondence between degrees and pixels
    assert(~any(xor(isnan(data.deg.vel),isnan(data.pix.vel))),'NaNs not same in data.deg.vel and data.pix.vel');
end

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
        blinkon (p) = sacon (qEnclosed);
        blinkoff(p) = sacoff(qEnclosed);
        
        % mark saccade as blink
        qIsBlink(qEnclosed) = true;
    else
        % all blinks should have been detected as saccades in the previous
        % steps
        warning('There is unphysiological movement in this trial that was not detected as a saccade. This might not have been dealt with properly although I''ll attempt. See sample %d, time %.1f ms\n',blinkon(p),blinkon(p)*1000/ETparams.samplingFreq);
        
        % see if eye position somehow plateaud during blink and thats
        % causing this trouble
        if nanmean(data.deg.Ele(blinkon(p):blinkoff(p)))> 10    % eye position more than 10 degrees off where it should be (again this is just for our horizontal tracking expt!!)
            % take previous and next saccade and take that as blink
            % interval
            previousSaccade = sacoff-blinkon(p);
            previousSaccadeidx = sum(previousSaccade<0);
            nextSaccade     = sacon-blinkoff(p);
            nextSaccadeidx  = sum(nextSaccade<0)+1;
            
            assert(previousSaccadeidx+1==nextSaccadeidx,'Couldn''t deal with potential blink...');
            
            blinkon (p) = sacon (previousSaccadeidx);
            blinkoff(p) = sacoff(    nextSaccadeidx);
            
            % mark saccades as blink
            qIsBlink([previousSaccadeidx nextSaccadeidx]) = true;
        else
            error('Couldn''t deal with potential blink...')
        end
    end
    
    if 0
        % debug, plot interval detected as blink
        figure(200),clf
        plot(data.deg.Ele(blinkon(p):blinkoff(p)))
        title(sprintf('blink %d',p))
        pause
    end
end

% multiple unphysiological segments might be enclosed by same saccade, in
% which case the same blink is detected multiple times. Remove duplicated
[~,idx]     = unique(blinkon);
blinkon     = blinkon (idx);
blinkoff    = blinkoff(idx);
% merge overlapping, this only happens when we dealt with trouble above....
% hack hacks need to be cleaned up
blink       = mergeIntervals(struct('on',blinkon,'off',blinkoff),[],0);
blinkon     = blink.on;
blinkoff    = blink.off;
% sanity checks
assert(all(blinkoff>blinkon) && all(blinkoff(1:end-1)<blinkon(2:end)))


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
    % -> We're dealing with a blink. Yes, thats rather arbitrary, but it
    %    did a good job for the data I was analyzing. You might want to do
    %    something else for another experiment (or another Eye tracker!).
    if nanmean(eye_pos)>0 && ...    % 1
            eye_max > 5 && ...      % 2
            eye_end < eye_max*.4    % 3
        
        % add info about blink
        blinkon  = [blinkon  sacon(p) ];
        blinkoff = [blinkoff sacoff(p)];
        % mark saccade as blink
        qIsBlink(p) = true;
        
        if 0
            % debug, plot interval detected as blink
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
    data.glissade   = replaceElementsInStruct(data.glissade,qRemoveGlissade,[]);
end

% then deal with the saccades
data.saccade    = replaceElementsInStruct(data.saccade,qIsBlink,[]);

% merge very close blinks.
data.blink      = mergeIntervals(data.blink, [], blinkMergeWindowSamples);

% replace with nan if wanted
if ETparams.blink.qReplaceWithNan
    % create boolean matrix given blink bounds
    qBlink = bounds2bool(data.blink.on+1,data.blink.off-1,length(data.deg.vel));    % remove one sample inwards as thats good for plotting and otherwise doesn't matter
    
    % remove data that is due to noise
    if ETparams.blink.qReplaceAllVelWithNan
        % do all velocity-based traces as well (we've found the relevant
        % fields above already!)
        data.deg    = replaceElementsInStruct(data.deg,qBlink,nan,velFieldsDeg);
        
        % if we have derivatives of eye position in pixels, throw the NaNs
        % in there as well
        if ETparams.data.qAlsoStoreandDiffPixels
            data.pix= replaceElementsInStruct(data.pix,qBlink,nan,velFieldsPix);
        end
    else
        % only the 2D velocity field 'vel'
        data.deg    = replaceElementsInStruct(data.deg,qBlink,nan,{'vel'});
        
        % if we have derivatives of eye position in pixels, throw the NaNs
        % in there as well
        if ETparams.data.qAlsoStoreandDiffPixels
            data.pix= replaceElementsInStruct(data.pix,qBlink,nan,{'vel'});
        end
    end
    
    % lastly, notify if more than 20% nan
    if sum(isnan(data.deg.vel))/length(data.deg.vel) > 0.20
        fprintf('Warning: This trial contains %.2f%% missing+blinks samples\n',sum(isnan(data.deg.vel))/length(data.deg.vel)*100);
        data.qNoiseTrial = true;
    else
        data.qNoiseTrial = false;
    end
end




% helper function for checking whether nans are in same positions for all
% traces
function velFields = checkNansSameForAllFields(data,datatype)

fn = fieldnames(data.(datatype));
% all position fields should have same NaNs, as should all velocity fields
% and all acceleration fields
qAccFields = cellfun(@(x) length(x)>2 && strcmp(x(1:3),'acc'),fn);
accFields  = fn(qAccFields);
for p=1:length(accFields)-1
    assert(~any(xor(isnan(data.(datatype).(accFields{p})),isnan(data.(datatype).(accFields{p+1})))),'NaNs not same in data.%1$s.%2$s and data.%1$s.%3$s',datatype,accFields{p},accFields{p+1});
end
qVelFields = cellfun(@(x) length(x)>2 && strcmp(x(1:3),'vel'),fn);
velFields  = fn(qVelFields);
for p=1:length(velFields)-1
    assert(~any(xor(isnan(data.(datatype).(velFields{p})),isnan(data.(datatype).(velFields{p+1})))),'NaNs not same in data.%1$s.%2$s and data.%1$s.%3$s',datatype,velFields{p},velFields{p+1});
end
posFields  = fn(~(qVelFields | qAccFields));
for p=1:length(posFields)-1
    assert(~any(xor(isnan(data.(datatype).(posFields{p})),isnan(data.(datatype).(posFields{p+1})))),'NaNs not same in data.%1$s.%2$s and data.%1$s.%3$s',datatype,posFields{p},posFields{p+1});
end