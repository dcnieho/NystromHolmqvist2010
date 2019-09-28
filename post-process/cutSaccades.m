function data = cutSaccades(data, ETparams,cutPosTraceMode,cutVelTraceMode,extraCut,skipFirstWindow)

% settings for the saccade cutting
% NB: depending on settings, the data traces velNoSac and posNoSac that are
% output to the data files might be the same as their unfiltered cousins..
% i.e., if you don't ask for any manipulation of these here.
%
% cutPosTraceMode: settings for the position trace:
% - 1: Use Ignace's method to reconstruct the pursuit that is hidden behind
%      the saccade, 2nd order interpolation of position.
% - 2: Reconstruct eye positions from desaccaded velocity signal (requires
%      cutVelTraceMode==1, see below)
% - 3: smooth out saccades in position domain by linear interpolation.
%
% cutVelTraceMode: settings for the position trace:
% - 1: Replace saccades by linearly interpolating velocity from saccade
%      onset to offset.
% - 2: Compute velocity trace from the reconstructed smooth pursuit signal
%      in position domain (cutPosTraceMode==1, see above).
%
% depending on your settings, some further inputs are needed:
% if cutPosTraceMode==1
%     % parameters for Ignace's method of removing pursuit from position
%     % trace
%     ETparams.Ignace.run     = 15;           % number of samples to use before and after saccade to fit polynomial to
%     ETparams.Ignace.offset  = 0;            % number of samples away from saccade start and end where RUN number of samples are taken.
% elseif cutPosTraceMode==2
%     % make plant for reconstructing position from velocity
%     Ts                      = 1/ETparams.samplingFreq;        % system sampling interval
%     sys                     = getSys(1,0,true,false,false);
%     ETparams.sysdt          = c2d(sys,Ts,'tustin');
% end
assert(cutVelTraceMode~=2 || cutPosTraceMode==1,'When reconstructing velocity from desaccaded eye position, desaccading must be run on the eye position trace')
assert(cutPosTraceMode~=2 || cutVelTraceMode==1,'When reconstructing position from desaccaded eye velocity, desaccading must be run on the eye velocity trace')
if cutPosTraceMode==1 && ~isfield(ETparams,'Ignace')
    % some default settings if user didn't provide them
    ETparams.Ignace.run     = 15;           % number of samples to use before and after saccade to fit polynomial to
    ETparams.Ignace.offset  = 0;
end
if cutPosTraceMode==2 && ~isfield(ETparams,'sysdt')
    % make a default plant for reconstructing position from velocity
    Ts              = 1/ETparams.samplingFreq;  % system sampling interval
    sys             = getSys(1,0,true,false,false);
    ETparams.sysdt  = c2d(sys,Ts,'tustin');
end

[data.deg.AziNoSac,...
 data.deg.EleNoSac,...
 data.deg.velNoSac,...
 data.deg.velAziNoSac,...
 data.deg.velEleNoSac] = cutSaccadesImplementation(data.deg.Azi,...
                                                  data.deg.Ele,...
                                                  data.deg.vel,...
                                                  data.deg.velAzi,...
                                                  data.deg.velEle,...
                                                  true,...
                                                  data.saccade,...
                                                  ETparams,cutPosTraceMode,cutVelTraceMode,extraCut,skipFirstWindow);

if isfield(data.pix,'vel')
    [data.pix.XNoSac,...
     data.pix.YNoSac,...
     data.pix.velNoSac,...
     data.pix.velXNoSac,...
     data.pix.velYNoSac] = cutSaccadesImplementation(data.pix.X,...
                                                    data.pix.Y,...
                                                    data.pix.vel,...
                                                    data.pix.velX,...
                                                    data.pix.velY,...
                                                    false,...
                                                    data.saccade,...
                                                    ETparams,cutPosTraceMode,cutVelTraceMode,extraCut,skipFirstWindow);
end



function [X,Y,vel,velX,velY] = cutSaccadesImplementation(X,Y,vel,velX,velY, qDeg,sac,ETparams,cutPosTraceMode,cutVelTraceMode,extraCut,skipFirstWindow)

% prepare parameters
% Number of samples at beginning of trial where we leave saccades
% untouched, if we reconstruct the position signal as well.
% Window length is in seconds
skipWindowSamples   = ceil(skipFirstWindow * ETparams.samplingFreq);

% stretch up the part cut around saccade onsets and offsets, if wanted
if nargin>2 && ~isempty(extraCut) && any(extraCut)
    sac.on  = sac.on  + ceil(extraCut(1)/1000 * ETparams.samplingFreq);
    sac.off = sac.off + ceil(extraCut(2)/1000 * ETparams.samplingFreq);
    
    % We could now have saccade starts before the end of the
    % previous saccade. prune/merge them here.
    sac = mergeIntervals(sac,[],0);
    
    % make sure first onset and last offset doesn't run out of the data
    if sac.on(1) < 1
        sac.on(1) = 1;
    end
    if sac.off(end) > length(X)
        sac.off(end) = length(X);
    end
end

sac.len= sac.off-sac.on+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with velocity
for p=1:length(sac.on) * (cutVelTraceMode==1)
    % per saccade, linearly interpolate eye velocity between saccade begin
    % point and end point
    if cutPosTraceMode==2 && sac.on(p) <= skipWindowSamples
        % special case: when reconstructing eye position, skip saccades
        % during the first second, they are probably to catch the blob as
        % it appears at the beginning of the trial. Otherwise we'll have a
        % constant offset in the data to deal with. We'll skip this part of
        % the data anyway for analysis, so no hurt incurred
        
        continue;
    end
    
    % replace with interpolated velocity
    [vel,velX,velY] = replaceIntervalVelocity(vel,velX,velY,Y,qDeg,sac.on(p),sac.off(p));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with position
switch cutPosTraceMode
    case 2
        % integrate velocities using plant
        if isnan(velX(1))
            sIdx = find(~isnan(velX),1);
        else
            sIdx = 1;
        end
        X(sIdx:end) = CanonicalDiscreteSSModel(ETparams.sysdt,velX(sIdx:end)).' + X(sIdx);
        Y(sIdx:end) = CanonicalDiscreteSSModel(ETparams.sysdt,velY(sIdx:end)).' + Y(sIdx);
    case 3
        % simply replace saccade with straight line between start and end
        % position. That isn't very useful if you want to remove the effect
        % of saccades from your data.
        for p=1:length(sac.on)
            X(sac.on(p):sac.off(p)) = linspace(X(sac.on(p)), X(sac.off(p)), sac.len(p));
            Y(sac.on(p):sac.off(p)) = linspace(Y(sac.on(p)), Y(sac.off(p)), sac.len(p));
        end
    case 1
        % Ignace interpolate - last as its long...
        for p=1:length(sac.on)
            if sac.on(p)<5 || sac.off(p)>length(X)-5
                % niet genoeg data om iets mee te kunnen doen, daar kunnen
                % we geen polynomial op fitten.
                % overslaan...
                continue;
            end
            
            % voorbeeld parameters en betekenis:
            % run   15: gebruik 15 samples voor en na saccade voor fitten,
            % offset 0: offset voor/na begin einde saccade waar de RUN lengte aan
            %           data gevonden wordt. Je kunt een beetje ruimte
            %           kiezen als je algoritme moeite heeft goede saccade
            %           beginnen en einden te vinden.
            soff    = ETparams.Ignace.offset+ETparams.Ignace.run;
            eoff    = ETparams.Ignace.offset;
            
            start1 = sac.on(p)-soff;
            if p==1
                % zorg dat we niet buiten de data grijpen
                start1 = max(start1,1);
            else
                % zorg dat we niet een stuk van een vorige saccade meenemen
                start1 = max(start1,sac.off(p-1));
            end
            eind1  = sac.on(p)-eoff;
            
            start2 = sac.off(p)+eoff;
            eind2  = sac.off(p)+soff;
            if p==length(sac.on)
                % zorg dat we niet buiten de data grijpen
                eind2 = min(eind2,length(X));
            else
                % zorg dat we niet een stuk van een volgende saccade
                % meenemen
                eind2 = max(eind2,sac.on(p+1));
            end
            
            X = interpSacIgnace(X,2,start1,eind1,start2,eind2,sac.on(p),sac.off(p)); % add ,true); for debug plots
            Y = interpSacIgnace(Y,2,start1,eind1,start2,eind2,sac.on(p),sac.off(p));
        
            
            % now, reconstruct velocity from desaccaded position, if wanted.
            % calc velocity during interpolated saccades and put that in
            if cutVelTraceMode==2
                on  = sac.on (p);
                off = sac.off(p);
                
                velX(on:off) = diff(X(on:off+1)) * ETparams.samplingFreq;   % TODO: change to using conv with [1 0 -1]
                velY(on:off) = diff(Y(on:off+1)) * ETparams.samplingFreq;
                
                vel(on:off)  = hypot(velX(on:off),velY(on:off));
            end
        end
end

if ismember(cutPosTraceMode,[1 2])
    % as most Ss show a preference to make saccades in a certain direction,
    % integrated position trace shows a trend. Remove trend by
    % least-squares procedure. Expect this trend to be approximately
    % opposite to trend of saccadic trace of course
    % See Collewijn & Tamminga - 1984 - Human smooth and saccadic eye
    % movements during voluntary pursuit of different target motions on
    % different backgrounds. J Physiol. 351, pp. 217-250.
    s = [1:length(X)].';
    t = regstats(X,s,'linear',{'beta'});
    X = X-s*t.beta(2);
    t = regstats(Y,s,'linear',{'beta'});
    Y = Y-s*t.beta(2);
end
