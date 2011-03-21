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

[data.saccade.on,data.saccade.off] = mergeSaccadesAndGlissades(data,ETparams,10);
sacon  = data.saccade.on;
sacoff = data.saccade.off;

% remove glissades, we fused them anyway
data = rmfield(data,'glissade');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% per saccade, linearly interpolate eye velocity between saccade begin
% point and end point
for p=1:length(sacon)
    % special case: skip saccades during the first second, they are 
    % probably to catch the blob as it appears at the beginning of the
    % trial. We skip this part of the data anyway for analysis
    if sacon(p) <= firstSecSamples
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
    [vel,velX,velY] = replaceSaccade(X,Y,vel,velX,velY,sacon(p),sacoff(p));
    
    if 0
        xc = CanonicalDiscreteSSModel(ETparams.sysdt,velX(sacon(p):sacoff(p))) + X(sacon(p));
        yc = CanonicalDiscreteSSModel(ETparams.sysdt,velY(sacon(p):sacoff(p))) + Y(sacon(p));
    end
end

% also interpolate NaNs that are still here
qNaN = isnan(velX);
if any(qNaN)
    fprintf('N NaN samples left: %d\n',sum(qNaN));
    [nanon,nanoff] = bool2bounds(qNaN);
    % gooi NaNnen aan begin en einde trial eruit, daar kunnen we niets mee
    if nanon(1)==1
        nanon(1)    = [];
        nanoff(1)   = [];
    end
    if ~isempty(nanoff) && nanoff(end)==length(velX)    % might be empty by now...
        nanon(end)  = [];
        nanoff(end) = [];
    end
    
    for p=1:length(nanon)
        % replace with interpolated velocity
        [vel,velX,velY] = replaceSaccade(X,Y,vel,velX,velY,nanon(p)-1,nanoff(p)+1); % pas indices aan, nanon(p) and nanoff(p) wijzen naar de eerste en laatste NaN in een serie
    end
    
    % show how many NaN we have left now, those cannot be handled
    fprintf(' -> N NaN samples left: %d\n',sum(isnan(velX)));
end



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
data.pix.Xfilt = CanonicalDiscreteSSModel(ETparams.sysdt,velX).' + X(1);
data.pix.Yfilt = CanonicalDiscreteSSModel(ETparams.sysdt,velY).' + Y(1);
data.pix.velfilt = vel;



%%% helpers
function [vel,velX,velY] = replaceSaccade(X,Y,vel,velX,velY,on,off)
% on and off are sample numbers (data indices) of saccade on- and offset
% respectively

% replace with interpolated velocity
vel(on:off)  = linspace(vel(on),vel(off),off-on+1);

% get saccade direction and use it to compute the X and Y velocity
% components
ang = atan2(Y(off)-Y(on),X(off)-X(on));
velX(on:off) = cos(ang)*vel(on:off);
velY(on:off) = sin(ang)*vel(on:off);