function data = detectAndRemoveBlinks(data,ETparams)
% Detects blinks and other un-physiological eye movements. They are
% replaced with nan or linear interpolation. Also detects some wobbly data
% as evidenced by noise in the pupil size data (almost always associated
% with noise in eye position, and especially velocity, data
% 
% Possible blinks are detect in one of two ways.
% 1. Based on peaks in change of pupil size velocity trace.
% 2. Velocity or acceleration trace shows speeds above what is
%    physiologically possible.


% prepare parameters
% if blink is followed by another blink by less than mergeWindow (in ms),
% they'll be merged
blinkMergeWindowSamples = ceil(ETparams.blink.mergeWindow           ./1000 * ETparams.samplingFreq);
minBlinkSamples         = ceil(ETparams.blink.minDur                ./1000 * ETparams.samplingFreq);
localNoiseWindowSamples = ceil(ETparams.blink.localNoiseWindowLength./1000 * ETparams.samplingFreq);


% First a sanity data integrity check: check NaN positions are in the same
% places in all traces, if not I must have introduced some bug somewhere.
% Not that its important for the correct functioning of this function by
% the way...
velFieldsDeg = checkNansSameForAllFields(data,'deg');

% check eye position in pixels as well (might only be position)
velFieldsPix = checkNansSameForAllFields(data,'pix');
% and check correspondence between degrees and pixels
assert(~any(xor(isnan(data.deg.Azi),isnan(data.pix.X))),'NaNs not same in data.deg.Azi and data.pix.X');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first get potential blinks, by samples above threshold and/or by
% detecting unphysiological eye movements

qBlink = false(size(data.deg.vel));
qNaN   = isnan(data.deg.vel);

% if pupil size is available, flag a blink when eye is closed/pupil lost
if isfield(data.pupil,'size')
    qBlink = qBlink | data.pupil.size==0;
end

% detect episodes above threshold in pupil size change trace
if bitget(uint8(ETparams.blink.detectMode),uint8(1))
    qBlink = qBlink | ...
             abs(data.pupil.dsize) > data.blink.peakDSizeThreshold;
end

% Detect possible blinks, episodes where the eyes move too fast too be
% physiologically possible
if bitget(uint8(ETparams.blink.detectMode),uint8(2))
    qBlink = qBlink | ...
             data.deg.vel > ETparams.blink.velocityThreshold |...
             data.deg.acc > ETparams.blink.accThreshold;
end

% find bounds of blinks (or tracker noise) as detected above, merge very
% close ones for easier processing, they'll merge anyway
[blink.on,blink.off]    = bool2bounds(qBlink);
blink                   = mergeIntervals(blink,[],blinkMergeWindowSamples,qNaN);

% if we have pupil size data, use it to refine blink onsets and offsets
if isfield(data.pupil,'dsize')
    % Process one possible blink at a time, refine the bounds
    % first, do refine like we do for saccades, if we're using peaks in the
    % pupil size change trace.
    dPSize  = abs(data.pupil.dsize);

    % make vector that will contain true where blinks are detected
    % this is used in some parts of the algorithm (e.g. offset velocity
    % threshold calculation) when it needs to analyze data that is not during
    % blinks or missing data.
    qBlinkOrNan   = qNaN;

    %%% Process one peak at the time.
    % Keep a counter here of how many peaks we have processed already. We need
    % to use a while loop instead of a for-loop as the length of the saccade
    % vector will change as we delete from it
    kk = 1;

    while kk <= length(blink.on)

        %----------------------------------------------------------------------  
        % Check the saccade peak samples to eliminate noise
        %----------------------------------------------------------------------

        % If the peak consists of <minPeakSamples consequtive samples, it it
        % probably noise, delete it. Except when pupil size reaches zero in
        % this short time, then its something we want to catch here
        if blink.off(kk)-blink.on(kk) < ETparams.blink.minPeakSamples && ...
           all(data.pupil.size(blink.on(kk):blink.off(kk))>0)
            blink.on (kk) = [];
            blink.off(kk) = [];
            continue;
        end

        % Detect blink start. Walk back from detected blink start to find
        % where the signal is below the onsetDSizeThreshold (mean+3*std)
        % and where the signal change is negative (which indicates a local
        % minimum in the signal).
        i = blink.on(kk);
        while i > 1 && ...                                          % make sure we don't run out of the data
              (isnan(dPSize(i)) || ...                              % and that we ignore nan data
               dPSize(i) > data.blink.onsetDSizeThreshold || ...    % keep searching until below onsetDSizeThreshold
               diff(dPSize(i-[0:1])) < 0 || ...                     % and signal change is negative
               any(data.pupil.size(i-[0 1])==0))                    % or pupil size is currently zero
            i = i-1;
        end
        blink.on(kk) = i;                                           % velocity minimum is last sample before "acceleration" sign change

        % Calculate local noise during period before the blink start (the
        % adaptive part), excluding data that is during known blink time or
        % missing.
        % This assumes the blink is preceded the same change in pupil size and
        % noise level as after the blink - its the best we can do.
        % starting from the already refined blink beginning we get
        % localNoiseWindowSamples samples before the blink, or as many as we
        % can get, that are not nan and not during blink.
        idx = find(~qBlinkOrNan(1:blink.on(kk)-1),localNoiseWindowSamples,'last');
        localVelNoise = dPSize(idx);
        localVelThresh= mean(localVelNoise) + 3*std(localVelNoise);

        % Check whether the local velocity threshold exceeds the peak velocity
        % threshold
        if ~isnan(localVelThresh) && localVelThresh < data.blink.peakDSizeThreshold
            data.blink.offsetDSizeThreshold(kk) = localVelThresh*0.3 + data.blink.onsetDSizeThreshold*0.7; % 30% local + 70% global
        else
            data.blink.offsetDSizeThreshold(kk) = data.blink.onsetDSizeThreshold;
        end

        % Detect blink end. Walk forward from detected blink end to find
        % where the signal is below the data.blink.offsetDSizeThreshold and
        % where the signal change is positive (which indicates a local minimum
        % in the signal)
        i = blink.off(kk);
        while i < length(dPSize) && ...                                 % make sure we don't run out of the data
              (isnan(dPSize(i)) || ...                                  % and that we ignore nan data
               dPSize(i) > data.blink.offsetDSizeThreshold(kk) || ...   % keep searching until below saccadeOffsetTreshold
               diff(dPSize(i+[0:1])) < 0 || ...                         % and signal change is positive
               any(data.pupil.size(i+[0 1])==0))                        % or pupil size is currently zero
            i = i+1;
        end
        blink.off(kk) = i;

        % now, delete all the next blinks that are enclosed by this potential
        % blink, i.e., delete any blink whose unrefined end is before
        % the end of this blink as they would converge to the same blink
        % interval. We do this now before the final checks below because if the
        % current blink is deleted by the below checks, any later blink
        % that will converge to the same interval would be deleted as well.
        while kk+1<=length(blink.off) &&...                            % make sure we don't run out of the data
              blink.off(kk+1) <= blink.off(kk)
            blink.on (kk+1) = [];
            blink.off(kk+1) = [];
            continue;
        end

        % Make sure the blink duration exceeds the minimum duration or
        % delete it. Except when pupil size reaches zero in this short
        % time, then its something we want to catch here
        if blink.off(kk)-blink.on(kk)+1 < minBlinkSamples&& ...
           all(data.pupil.size(blink.on(kk):blink.off(kk))>0)
            blink.on (kk) = [];
            blink.off(kk) = [];
            continue;
        end

        %%%%
        % Done. All the above criteria are fulfilled, we've got a blink.
        %%%%
        % flag it in the trace
        qBlinkOrNan(blink.on(kk):blink.off(kk)) = true;

        % increase counter, process next peak
        kk = kk+1;
    end
end

% Then, attempt to further refine blinks by noting that many of these
% unphysiological eye movements have already been detected as saccades
% (they usually have very large velocity!). So use detected saccade on- and
% offsets as blink on- and offsets if the first saccade onset is earlier or
% the last offset later (crap in velocity trace is sometimes a bit more
% spread out)
sacon  = data.saccade.on;
sacoff = data.saccade.off;
% correct sacoff to glissade end time, if any. we wnat the blink to extend
% all the way to glissade offset
if isfield(data,'glissade')
    qHaveGliss = ismember(data.saccade.off,data.glissade.on);
    sacoff(qHaveGliss) = data.glissade.off;
end
for p = 1:length(blink.on)
    % any saccade on/off during blink interval, or is blink interval enclosed by saccade on off?
    qSac = (sacon>=blink.on(p) & sacon<=blink.off(p)) | (sacoff>=blink.on(p) & sacoff<=blink.off(p)) | (blink.on(p)>sacon & blink.on(p)<sacoff);
    
    blink.on (p) = min(blink.on(p) , minOrInf (sacon (qSac)));
    blink.off(p) = max(blink.off(p), maxOrNInf(sacoff(qSac)));
end

% multiple unphysiological segments might be enclosed by same saccade, or
% otherwise get refined to same onsets/offsets, in which case the same
% blink is detected multiple times. Remove duplicated
[~,idx]     = unique(blink.on);
blink.on    = blink.on (idx);
blink.off   = blink.off(idx);

% merge very close blinks.
blink       = mergeIntervals(blink, [], blinkMergeWindowSamples,qNaN);

% Make sure pupil size trace shows a peak. If fast unidirectional
% change, we caught something else accidentally, such as a temporarily
% drooping eyelid, or an actual sudden change in pupil size. Check
% proportion of positive change in pupil size as a function of total
% non-zero changes in pupil size. if close to 0 or close to 1, skip
% this
if isfield(data.pupil,'dsize')
    for p = length(blink.on):-1:1
        propPos = sum(data.pupil.dsize(blink.on(p):blink.off(p))>0) / sum(data.pupil.dsize(blink.on(p):blink.off(p))~=0);
        if propPos<.05 || propPos>.95
            blink.on (p) = [];
            blink.off(p) = [];
            continue;
        end
    end
end

% sanity check
assert(all(blink.off>blink.on) && all(blink.off(1:end-1)<blink.on(2:end)))

% build information about blinks
[blink.on,idx]  = sort(blink.on);
data.blink.on   = blink.on;
data.blink.off  = blink.off(idx);

% now remove saccades that were flagged as blink
% any saccade on/off during blink interval, or is blink interval enclosed
% by saccade on/off? Off is glissade off if the saccade is followed by a
% glissade, see above. This is what we want.
qIsBlink = false(size(sacon));
for p=1:length(blink.on)
    qIsBlink = qIsBlink | (sacon>=blink.on(p) & sacon<=blink.off(p)) | (sacoff>=blink.on(p) & sacoff<=blink.off(p)) | (blink.on(p)>sacon & blink.on(p)<sacoff);
end

% first remove their corresponding glissades, if any
if isfield(data,'glissade')
    qRemoveGlissade = ismember(data.glissade.on,data.saccade.off(qIsBlink));
    data.glissade   = replaceElementsInStruct(data.glissade,qRemoveGlissade,[]);
end

% then deal with the saccades
data.saccade    = replaceElementsInStruct(data.saccade,qIsBlink,[]);

% replace by linear interpolation or with nan if wanted
nNaN = sum(isnan(data.deg.vel));
nBlink = sum(blink.off-blink.on+1);
if ETparams.blink.qReplaceWithInterp
    % adjust indices, blink.on(p) and blink.off(p) point to first and last
    % samples of a blink
    blon  = blink.on-1;  blon (blon <1                   ) = 1;
    bloff = blink.off+1; bloff(bloff>length(data.deg.Azi)) = length(data.deg.Azi);
    bllen = bloff-blon+1;
    
    for p=1:length(blon)
        % position
        data.deg.Azi(blon(p):bloff(p)) = linspace(data.deg.Azi(blon(p)), data.deg.Azi(bloff(p)), bllen(p));
        data.deg.Ele(blon(p):bloff(p)) = linspace(data.deg.Ele(blon(p)), data.deg.Ele(bloff(p)), bllen(p));
        data.pix.X  (blon(p):bloff(p)) = linspace(data.pix.X  (blon(p)), data.pix.X  (bloff(p)), bllen(p));
        data.pix.Y  (blon(p):bloff(p)) = linspace(data.pix.Y  (blon(p)), data.pix.Y  (bloff(p)), bllen(p));
    
        % velocity
        [data.deg.vel,data.deg.velAzi,data.deg.velEle] = ...
            replaceIntervalVelocity(data.deg.vel,data.deg.velAzi,data.deg.velEle,data.deg.Ele,...
                                    true,blon(p),bloff(p));
        if isfield(data.pix,'velX')
            [data.pix.vel,data.pix.velX,data.pix.velY] = ...
                replaceIntervalVelocity(data.pix.vel,data.pix.velX,data.pix.velY,[],...
                                        false,blon(p),bloff(p));
        elseif isfield(data.pix,'vel')
            data.pix.vel(blon(p):bloff(p)) = linspace(data.pix.vel(blon(p)), data.pix.vel(bloff(p)), bllen(p));
        end

        % acceleration (abuse replaceintervalVel, keep qDeg input always fals so simple hypot() is used for 2D acc)
        [data.deg.acc,data.deg.accAzi,data.deg.accEle] = ...
            replaceIntervalVelocity(data.deg.acc,data.deg.accAzi,data.deg.accEle,[],...
                                    false,blon(p),bloff(p));
        if isfield(data.pix,'accX')
            [data.pix.acc,data.pix.accX,data.pix.accY] = ...
                replaceIntervalVelocity(data.pix.acc,data.pix.accX,data.pix.accY,[],...
                                        false,blon(p),bloff(p));
        elseif isfield(data.pix,'acc')
            data.pix.acc(blon(p):bloff(p)) = linspace(data.pix.acc(blon(p)), data.pix.acc(bloff(p)), bllen(p));
        end
    end
    
    % lastly, notify how much blinks
    if nBlink>0
        fprintf('  N blink samples: %d (%.2f%%)\n',nBlink,nBlink/length(data.deg.vel)*100);
    end
    
elseif ETparams.blink.qReplaceVelWithNan
    % create boolean matrix given blink bounds
    qBlink = bounds2bool(data.blink.on+1,data.blink.off-1,length(data.deg.vel));    % remove one sample inwards as thats good for plotting and otherwise doesn't matter
    
    % remove data that is due to noise
    % do all velocity-based traces as well (we've found the relevant
    % fields above already!)
    data.deg    = replaceElementsInStruct(data.deg,qBlink,nan,velFieldsDeg);
    
    % if we have derivatives of eye position in pixels, throw the NaNs
    % in there as well
    if ETparams.data.qAlsoStoreandDiffPixels
        data.pix= replaceElementsInStruct(data.pix,qBlink,nan,velFieldsPix);
    end
end
    
% lastly, notify if more than 20% nan
if (nBlink+nNaN)/length(data.deg.vel) > 0.20
    fprintf('Warning: This trial contains %.2f%% missing+blinks samples\n',(nBlink+nNaN)/length(data.deg.vel)*100);
    data.qNoiseTrial = true;
else
    data.qNoiseTrial = false;
end




% helper function for checking whether nans are in same positions for all
% traces
function [velFields,posFields,accFields] = checkNansSameForAllFields(data,datatype)

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

function out = minOrInf(in)
if isempty(in)
    out = inf;
else
    out = min(in(:));
end

function out = maxOrNInf(in)
if isempty(in)
    out = -inf;
else
    out = max(in(:));
end