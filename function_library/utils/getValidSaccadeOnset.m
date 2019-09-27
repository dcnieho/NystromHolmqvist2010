function varargout = getValidSaccadeOnset(dat,eyePosSamps,displacementThresh,eyePosThresh,maxSacDur,normalFields,noGapFields,target)

% skip saccades that are:
% 1. position peaks (so no actual displacement during saccade): more
% generally, saccades need to have a minimum horizontal (or 2D)
% displacement to count
% 2. eye position before saccade (not at start trial) needs to be within
% threshold from center of screen
% 3. "saccade" not too long, must be some noise crap
% 4. saccade that has missing directly after its onset, or directly before
% its offset
%
% 1 means we should skip the event in question and look at next
% 2 and 3 mean invalid trial: return NaN. For 4, return nan if missign
% directly after onset, return nan only for duration if missing directionly
% before offset

qTarget = nargin>=8;

assert(nargout==length(normalFields)+length(noGapFields)+2*qTarget,'number of fields given does not match number of requested output arguments')

% no data for this trial
if isempty(dat)
    [varargout{1:nargout}] = deal(nan);
    return;
end

saccade = dat.saccade;
azi     = dat.deg.Azi;
ele     = dat.deg.Ele;
time    = dat.time;

% see if any saccades in the first place
if isempty(saccade.on)
    [varargout{1:nargout}] = deal(nan);
    return;
end

% start at first saccade after target onset
start = find(time>=0,1);
sacIdx = find(dat.saccade.on>start,1);
if isempty(sacIdx)
    [varargout{1:nargout}] = deal(nan);
    return;
end

% check 1: position peaks. get pos before and after saccade (mean of
% eyePosSamps before and after saccade) and check if displacement exceeds
% threshold
while true
    if displacementThresh(1)==1
        % check horizontal
        sidx            = [-eyePosSamps:0]+saccade.on (sacIdx); sidx(sidx<1) = [];
        startPointAzi   = mean(azi(sidx),'omitnan');
        eidx            = [0:eyePosSamps] +saccade.off(sacIdx); eidx(eidx>length(azi)) = [];
        endPointAzi     = mean(azi(eidx),'omitnan');
        amp             = abs(endPointAzi-startPointAzi);
    elseif displacementThresh(1)==2
        % check vertical
        error('not implemented');
    elseif displacementThresh(1)==3
        % check 2D amplitude
        amp = saccade.amplitude(sacIdx);
    end
    if amp < displacementThresh(2)
        sacIdx = sacIdx+1;
        if sacIdx>length(saccade.on)
            [varargout{1:nargout}] = deal(nan);
            return;
        end
    else
        break;
    end
end
    
% check 2: eye position before saccade close enough to center of screen
sidx = [-eyePosSamps:0]+saccade.on (sacIdx); sidx(sidx<1) = [];
startPointAzi = mean(azi(sidx),'omitnan');
startPointEle = mean(ele(sidx),'omitnan');
if calcAmplitudeFick(0,0, startPointAzi,startPointEle)>eyePosThresh
    [varargout{1:nargout}] = deal(nan);
    return;
end

% check 3: saccade not too long
if saccade.duration(sacIdx) > maxSacDur
    [varargout{1:nargout}] = deal(nan);
    return;
end
    
% check 4: onset and offset are not at the edge of missing velocity data
qNoDuration = false;
if isnan(dat.deg.vel(saccade.on(sacIdx)+1))
    % useless saccade as onset is not reliable
    [varargout{1:nargout}] = deal(nan);
    return;
elseif isnan(dat.deg.vel(saccade.off(sacIdx)-1))
    % saccade duration useless as offset not reliable
    qNoDuration = true;
end



% get requested info for selected saccade
varargout = cell(1,nargout);
varargout(1:length(normalFields)) = cellfun(@(x) saccade.(x)(sacIdx),normalFields,'uni',false);

% if no missing during saccade (even single points), then also get the no
% gap fields, else return nan for those
if any(isnan(azi(saccade.on(sacIdx):saccade.off(sacIdx))))
    [varargout{length(normalFields)+[1:length(noGapFields)]}] = deal(nan);
else
     varargout(length(normalFields)+[1:length(noGapFields)]) = cellfun(@(x) saccade.(x)(sacIdx),noGapFields,'uni',false);
end

if qNoDuration
    qDur = strcmp([normalFields noGapFields],'duration');
    [varargout{qDur}] = deal(nan);
end

if qTarget
    % see if we have info about where target is. then we can see if second
    % saccade is a corrective saccade that bring eye closer to target
    % output size of error reduction and onset t
    if sacIdx<length(saccade.on)
        % get saccade onset and offset position
        sacIdx          = sacIdx+1;
        sidx            = [-eyePosSamps:0]+saccade.on (sacIdx); sidx(sidx<1) = [];
        startPointAzi   = mean(azi(sidx),'omitnan');
        startPointEle   = mean(ele(sidx),'omitnan');
        eidx            = [0:eyePosSamps] +saccade.off(sacIdx); eidx(eidx>length(azi)) = [];
        endPointAzi     = mean(azi(eidx),'omitnan');
        endPointEle     = mean(ele(eidx),'omitnan');
        
        % get target position
        target = target(dat.file.trial+1,4:5);
        % calc error reduction
        errReduc    = calcAmplitudeFick(target(1),target(2), startPointAzi,startPointEle) - calcAmplitudeFick(target(1),target(2), endPointAzi,endPointEle);
        if errReduc>0
            varargout{end-1} = errReduc;
            varargout{end}   = saccade.onPrecise(sacIdx);
        else
            [varargout{end-1:end}] = deal(nan);
        end
    else
        [varargout{end-1:end}] = deal(nan);
    end
end

% last, correct onset and offset fields for t==0 point, convert to time
outFields = [normalFields noGapFields];
toCorrect = find(ismember(outFields,{'on','off','onPrecise'}));
for p=1:length(toCorrect)
    varargout{toCorrect(p)} = interp1(1:length(time),time,varargout{toCorrect(p)});
end
% also correct duration
toCorrect = find(ismember(outFields,{'duration'}));
for p=1:length(toCorrect)
    varargout{toCorrect(p)} = interp1(1:length(time),time,saccade.off(sacIdx))-interp1(1:length(time),time,saccade.on(sacIdx));
end
