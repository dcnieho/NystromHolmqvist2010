function data = cutSaccades(data,datatype,ETparams, qReconstructPos,skipFirstWindow)

% prepare parameters
% Number of samples at beginning of trial where we leave saccades
% untouched, if we reconstruct the position signal as well.
% Window length is in seconds
skipWindowSamples   = ceil(skipFirstWindow * ETparams.samplingFreq);

% get eye positions in pixels/degree
if strcmp(datatype,'pix')
    X   = data.pix.X;
    Y   = data.pix.Y;
elseif strcmp(datatype,'deg')
    X   = data.deg.Azi;
    Y   = data.deg.Ele;
end

% get eye velocities in pixels/degree
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
% per saccade, linearly interpolate eye velocity between saccade begin
% point and end point
for p=1:length(sacon)
    % special case: skip saccades during the first second, they are
    % probably to catch the blob as it appears at the beginning of the
    % trial. We skip this part of the data anyway for analysis
    % DN: this is only important for reconstructing eye position, so only
    % use it if we want that
    if qReconstructPos && sacon(p) <= skipWindowSamples
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
        [vel,velX,velY] = replaceIntervalVelocity(vel,velX,velY,sacon(p),sacoff(p));
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