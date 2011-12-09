function data = cutSaccades(data, ETparams,cutPostrace,skipFirstWindow)

[data.deg.AziFilt,...
 data.deg.EleFilt,...
 data.deg.velFilt,...
 data.deg.velAziFilt,...
 data.deg.velEleFilt] = cutSaccadesImplementation(data.deg.Azi,...
                                                  data.deg.Ele,...
                                                  data.deg.vel,...
                                                  data.deg.velAzi,...
                                                  data.deg.velEle,...
                                                  data.saccade,...
                                                  ETparams,cutPostrace,skipFirstWindow);

if isfield(data.pix,'vel')
    [data.pix.XFilt,...
     data.pix.YFilt,...
     data.pix.velFilt,...
     data.pix.velXFilt,...
     data.pix.velYFilt] = cutSaccadesImplementation(data.pix.X,...
                                                    data.pix.Y,...
                                                    data.pix.vel,...
                                                    data.pix.velX,...
                                                    data.pix.velY,...
                                                    data.saccade,...
                                                    ETparams,cutPostrace,skipFirstWindow);
end



function [X,Y,vel,velX,velY] = cutSaccadesImplementation(X,Y,vel,velX,velY, sac,ETparams,cutPostrace,skipFirstWindow)

% prepare parameters
% Number of samples at beginning of trial where we leave saccades
% untouched, if we reconstruct the position signal as well.
% Window length is in seconds
skipWindowSamples   = ceil(skipFirstWindow * ETparams.samplingFreq);

sac.len= sac.off-sac.on+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% per saccade, linearly interpolate eye velocity between saccade begin
% point and end point
for p=1:length(sac.on)
    % special case: skip saccades during the first second, they are
    % probably to catch the blob as it appears at the beginning of the
    % trial. We skip this part of the data anyway for analysis
    % DN: this is only important for reconstructing eye position, so only
    % use it if we want that
    if cutPostrace==1 && sac.on(p) <= skipWindowSamples
        if any(isnan(vel(sac.on(p):sac.off(p))))
            % if there is some NaN during this first interval, create a
            % position that linearly interpolates between begin and end
            tempX = linspace(X(sac.on(p)), X(sac.off(p)), sac.len(p));
            tempY = linspace(Y(sac.on(p)), Y(sac.off(p)), sac.len(p));
            
            % replace velocity signal with this, creating the fake saccade
            % we don't have to care this messes up the spectrum a bit as we
            % skip the first 5 seconds anyway
            velX(sac.on(p):sac.off(p)) = diff(tempX(1:2))*ETparams.samplingFreq;    % constant velocity, so just compute it once
            velY(sac.on(p):sac.off(p)) = diff(tempY(1:2))*ETparams.samplingFreq;
            vel(sac.on(p):sac.off(p)) = hypot(diff(tempX(1:2)), diff(tempY(1:2)))*ETparams.samplingFreq;
        end
    else
        
        % replace with interpolated velocity
        [vel,velX,velY] = replaceIntervalVelocity(vel,velX,velY,sac.on(p),sac.off(p));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with position
switch cutPostrace
    case 1
        % integrate velocities using plant
        X = CanonicalDiscreteSSModel(ETparams.sysdt,velX).' + X(1);
        Y = CanonicalDiscreteSSModel(ETparams.sysdt,velY).' + Y(1);
    case 2
        % simply linearly interpolate position
        for p=1:length(sac.on)
            X(sac.on(p):sac.off(p)) = linspace(X(sac.on(p)), X(sac.off(p)), sac.len(p));
            Y(sac.on(p):sac.off(p)) = linspace(Y(sac.on(p)), Y(sac.off(p)), sac.len(p));
        end
end