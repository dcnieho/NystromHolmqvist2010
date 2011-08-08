function data = cutSaccades(data,ETparams,qReconstructPos)

% prepare parameters
firstSecSamples = 1 * ETparams.samplingFreq;    % number of samples in one second (looks stupid, but if we later decide to want a different interval, we can change things here...)

% get eye positions in pixels
X       = data.pix.X;
Y       = data.pix.Y;

% get eye velocities in pixels
vel     = data.pix.vel;
velX    = -data.pix.velX;   % not sure why, but the Savitzky-Golay filter gives me the wrong sign for the component velocities
velY    = -data.pix.velY;
%velX(isnan(velX)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocess on and offsets:
% - fuse glissade after saccade
% - fuse two saccades with little time between them
% - the two above also incorporates: fuse saccade that starts at glissade offset

[data.saccade.on,data.saccade.off] = mergeSaccadesAndGlissades(data,ETparams,10);
sacon  = data.saccade.on;
sacoff = data.saccade.off;

% remove glissades, we fused them anyway
data = rmfield(data,'glissade');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step one: we want to deal with all the nan in the data, also the position
% data
qNaN = isnan(vel);
if any(qNaN)
    fprintf('N NaN samples: %d\n',sum(qNaN));
    [nanon,nanoff] = bool2bounds(qNaN);
    % gooi NaNnen aan begin en einde trial eruit, daar kunnen we niets mee
    if nanon(1)==1
        nanon(1)    = [];
        nanoff(1)   = [];
    end
    if ~isempty(nanoff) && nanoff(end)==length(vel)     % might be empty by now...
        nanon(end)  = [];
        nanoff(end) = [];
    end
    
    for p=1:length(nanon)
        qDuringSac = nanon(p)>=sacon & nanoff(p)<=sacoff;
        if any(qDuringSac)
            % if nan is during saccade, use those as start and end points
            % (probably blink)
            assert(sum(qDuringSac)==1)  % anything else would be ridiculous!
            on  = sacon (qDuringSac);
            off = sacoff(qDuringSac);
            sacon (qDuringSac) = [];
            sacoff(qDuringSac) = [];
        else
            % pas indices aan, nanon(p) and nanoff(p) wijzen naar de eerste
            % en laatste NaN in een serie
            on  = nanon(p)-1;
            off = nanoff(p)+1;
        end
        % replace with interpolated velocity
        if qReconstructPos
            [vel,velX,velY] = replaceSaccade(X,Y,vel,velX,velY,on,off);
        else
            vel = replaceSaccade(X,Y,vel,velX,velY,on,off);
        end
    end
    
    % replace original vel with this one as we'll need one with nans
    % removed
    data.pix.vel = vel;
    
    % show how many NaN we have left now, those cannot be handled
    fprintf(' -> N NaN samples left: %d\n',sum(isnan(vel)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% per saccade, linearly interpolate eye velocity between saccade begin
% point and end point
for p=1:length(sacon)
    % special case: skip saccades during the first second, they are
    % probably to catch the blob as it appears at the beginning of the
    % trial. We skip this part of the data anyway for analysis
    if 0 && sacon(p) <= firstSecSamples
        if any(isnan(vel(sacon(p):sacoff(p))))
            % if there is some NaN during this first interval, create a
            % position that linearly interpolates between begin and end
            tempX = linspace(X(sacon(p)), X(sacoff(p)), sacoff(p)-sacon(p)+1);
            tempY = linspace(Y(sacon(p)), Y(sacoff(p)), sacoff(p)-sacon(p)+1);
            
            % replace velocity signal with this, creating the fake saccade
            % we don't have to care this messes up the spectrum a bit as we
            % skip the first 5 seconds anyway
            velX(sacon(p):sacoff(p)) = diff(tempX(1:2))*ETparams.samplingFreq;    % constant velocity, so just compute it once
            velY(sacon(p):sacoff(p)) = diff(tempY(1:2))*ETparams.samplingFreq;
            vel(sacon(p):sacoff(p)) = hypot(diff(tempX(1:2)), diff(tempY(1:2)))*ETparams.samplingFreq;
        end
        continue;
    end
    
    % replace with interpolated velocity
    if qReconstructPos
        [vel,velX,velY] = replaceSaccade(X,Y,vel,velX,velY,sacon(p),sacoff(p));
    else
        vel = replaceSaccade(X,Y,vel,velX,velY,sacon(p),sacoff(p));
    end
end

if qReconstructPos
    % plant version
    data.pix.Xfilt = CanonicalDiscreteSSModel(ETparams.sysdt,velX).' + X(1);
    data.pix.Yfilt = CanonicalDiscreteSSModel(ETparams.sysdt,velY).' + Y(1);
end
data.pix.velfilt = vel;



%%% helpers
function [vel,velX,velY] = replaceSaccade(X,Y,vel,velX,velY,on,off)
% on and off are sample numbers (data indices) of saccade on- and offset
% respectively

% replace with interpolated velocity
vel(on:off)  = linspace(vel(on),vel(off),off-on+1);

if nargout>1
    % get saccade direction and use it to compute the X and Y velocity
    % components
    ang = atan2(Y(off)-Y(on),X(off)-X(on));
    velX(on:off) = cos(ang)*vel(on:off);
    velY(on:off) = sin(ang)*vel(on:off);
end