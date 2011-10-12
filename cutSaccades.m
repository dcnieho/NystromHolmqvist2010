function data = cutSaccades(data,datatype,ETparams,qReconstructPos)

% prepare parameters
firstSecSamples = 1 * ETparams.samplingFreq;    % number of samples in one second (looks stupid, but if we later decide to want a different interval, we can change things here...)

% get eye positions in pixels/degree
if strcmp(datatype,'pix')
    X   = data.pix.X;
    Y   = data.pix.Y;
elseif strcmp(datatype,'deg')
    X   = data.deg.Azi;
    Y   = data.deg.Ele;
end

% get eye velocities in pixels
vel     = data.(datatype).vel;
if strcmp(datatype,'pix')
    velX    = data.pix.velX;
    velY    = data.pix.velY;
elseif strcmp(datatype,'deg')
    velX    = data.deg.velAzi;
    velY    = data.deg.velEle;
end

sacon  = data.saccade.on;
sacoff = data.saccade.off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step one: we want to deal with all the nan in the data.
% This is getting rid of blinks and such...
qNaN = isnan(vel);
if any(qNaN)
    fprintf('N NaN samples: %d\n',sum(qNaN));
    [nanon,nanoff] = bool2bounds(qNaN);
    % gooi NaNnen gevonden aan begin en einde trial eruit, daar kunnen we
    % niets mee...
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
        [vel,velX,velY] = replaceSaccade(vel,velX,velY,on,off);
    end
    
    % show how many NaN we have left now, those cannot be handled
    fprintf(' -> N NaN samples left: %d\n',sum(isnan(vel)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % replace original vel with this one as we'll need one with nans
    % removed
    data.(datatype).vel = vel;
    if strcmp(datatype,'pix')
        data.pix.velX = velX;
        data.pix.velY = velY;
    elseif strcmp(datatype,'deg')
        data.deg.velAzi = velX;
        data.deg.velEle = velY;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% per saccade, linearly interpolate eye velocity between saccade begin
% point and end point
for p=1:length(sacon)
    % special case: skip saccades during the first second, they are
    % probably to catch the blob as it appears at the beginning of the
    % trial. We skip this part of the data anyway for analysis
    % DN: this is only important for reconstructing eye position, so only
    % use it if we want that
    if qReconstructPos && sacon(p) <= firstSecSamples
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
    else
        
        % replace with interpolated velocity
        [vel,velX,velY] = replaceSaccade(vel,velX,velY,sacon(p),sacoff(p));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% store filtered data
data.(datatype).velFilt = vel;
if strcmp(datatype,'pix')
    data.pix.velXFilt = velX;
    data.pix.velYFilt = velY;
elseif strcmp(datatype,'deg')
    data.deg.velAziFilt = velX;
    data.deg.velEleFilt = velY;
end

if qReconstructPos
    % plant version
    XFilt = CanonicalDiscreteSSModel(ETparams.sysdt,velX).' + X(1);
    YFilt = CanonicalDiscreteSSModel(ETparams.sysdt,velY).' + Y(1);
    
    if strcmp(datatype,'pix')
        data.pix.XFilt = XFilt;
        data.pix.YFilt = YFilt;
    elseif strcmp(datatype,'deg')
        data.deg.AziFilt = XFilt;
        data.deg.EleFilt = YFilt;
    end
end



%%% helpers
function [vel,velX,velY] = replaceSaccade(vel,velX,velY,on,off)
% on and off are sample numbers (data indices) of saccade on- and offset
% respectively

npoint = off-on+1;

% components: replace with linearly interpolated velocity
velX(on:off) = linspace(velX(on),velX(off),npoint);
velY(on:off) = linspace(velY(on),velY(off),npoint);

% calculate interpolate 2D velocity. In effect this is now interpolated
% with a bicubic spline. Thats fine, good even as no edges are introduced
% into the data, as we care most about the component velocities in the
% situations I can think of.
vel(on:off)  = hypot(velX(on:off),velY(on:off));