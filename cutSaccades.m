function data = cutSaccades(data,ETparams)

% prepare parameters
firstSecSamples = 1 * ETparams.samplingFreq;    % number of samples in one second (looks stupid, but if we later decide to want a different interval, we can change things here...)

% get eye positions in pixels
X       = data.pix.X;
Y       = data.pix.Y;

% get eye velocities in pixels
vel     = data.pix.vel;
velX    = -data.pix.velX;   % not sure why, but the Savitzky-Golay filter gives me th wrong sign for the component velocities
velY    = -data.pix.velY;
%velX(isnan(velX)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocess on and offsets:
% - fuse glissade after saccade
% - fuse two saccades with little time between them
% - the two above also incorporates: fuse saccade that starts at glissade offset

[data.saccade.on,data.saccade.off] = MergeSaccadesAndGlissades(data,ETparams,10);
sacon  = data.saccade.on;
sacoff = data.saccade.off;

% remove glissades, we fused them anyway
data = rmfield(data,'glissade');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% per saccade, linearly interpolate eye velocity between saccade begin
% point and end point
for p=1:length(sacon)
    % special case: skip first saccade if it is during the first second, is
    % probably to catch the blob as it appears at the beginning of the
    % trial
    if p==1 && sacon(p) <= firstSecSamples
        continue;
    end
    
    % replace with interpolated velocity
    vel(sacon(p):sacoff(p))  = linspace(vel(sacon(p)),vel(sacoff(p)),sacoff(p)-sacon(p)+1);
    
    % get saccade direction and use it to compute the X and Y velocity
    % components
    ang = atan2(Y(sacoff(p))-Y(sacon(p)),X(sacoff(p))-X(sacon(p)));
    velX(sacon(p):sacoff(p)) = cos(ang)*vel(sacon(p):sacoff(p));
    velY(sacon(p):sacoff(p)) = sin(ang)*vel(sacon(p):sacoff(p));
    
    if 0
        xc = CanonicalDiscreteSSModel(ETparams.sysdt,velX(sacon(p):sacoff(p))) + X(sacon(p));
        yc = CanonicalDiscreteSSModel(ETparams.sysdt,velY(sacon(p):sacoff(p))) + Y(sacon(p));
    end
end

% TODO, also interpolate NaNs that are still here
fprintf('N NaN samples left: %d\n',sum(isnan(velX)));

if 0
% reconstruct eye position from the altered velocity signal
% x_{k+1} = x_k + v_k*delta_t
int = 1/ETparams.samplingFreq;
for p=1:length(velX)
    if p==1
        ddxx(p) = X(1);
    else
        ddxx(p) = ddxx(p-1)+velX(p-1)*int;
    end
end
end

% plant version
data.pix.Xfilt = CanonicalDiscreteSSModel(ETparams.sysdt,velX) + X(1);
data.pix.Yfilt = CanonicalDiscreteSSModel(ETparams.sysdt,velY) + Y(1);
data.pix.velfilt = vel;