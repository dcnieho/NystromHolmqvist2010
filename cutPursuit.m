function data = cutPursuit(data, ETparams,skipFirstWindow)

[data.deg.AziSac,...
 data.deg.EleSac,...
 data.deg.velSac,...
 data.deg.velAziSac,...
 data.deg.velEleSac] = cutPursuitImplementation(data.deg.Azi,...
                                                  data.deg.Ele,...
                                                  data.deg.vel,...
                                                  data.deg.velAzi,...
                                                  data.deg.velEle,...
                                                  true,...
                                                  data.saccade,...
                                                  ETparams,skipFirstWindow);

if isfield(data.pix,'vel')
    [data.pix.XSac,...
     data.pix.YSac,...
     data.pix.velSac,...
     data.pix.velXSac,...
     data.pix.velYSac] = cutPursuitImplementation(data.pix.X,...
                                                    data.pix.Y,...
                                                    data.pix.vel,...
                                                    data.pix.velX,...
                                                    data.pix.velY,...
                                                    false,...
                                                    data.saccade,...
                                                    ETparams,skipFirstWindow);
end



function [X,Y,vel,velX,velY] = cutPursuitImplementation(X,Y,vel,velX,velY, qDeg,sac,ETparams,skipFirstWindow)

% prepare parameters
% Number of samples at beginning of trial where we leave data
% untouched, if we reconstruct the position signal as well.
% Window length is in seconds
skipWindowSamples   = ceil(skipFirstWindow * ETparams.samplingFreq);

% get parts of data that is not saccades (and thus counts as pursuit)
nSamp  = length(X);
[pur.on,pur.off] = invertBounds(sac.on,sac.off,nSamp);
pur.len          = pur.off-pur.on+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with velocity
for p=1:length(pur.on)
    % per saccade, linearly interpolate eye velocity between saccade begin
    % point and end point
    if pur.on(p) <= skipWindowSamples
        % skip data during the first second, they are probably to catch the
        % blob as it appears at the beginning of the trial. Otherwise we'll
        % have a constant offest in the data to deal with. We'll skip this
        % part of the data anyway for analysis, so no hurt incurred
        
        continue;
    end
    
    % replace with 0 velocity
    [vel,velX,velY] = replaceIntervalVelocity(vel,velX,velY,Y,qDeg,pur.on(p),pur.off(p),0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with position
for p=1:length(pur.on)
    % simply hold value of sample before pursuit onset. Then shift trace
    % after to match that value
    if pur.on(p) <= skipWindowSamples
        % skip data during the first second, they are probably to catch the
        % blob as it appears at the beginning of the trial. Otherwise we'll
        % have a constant offest in the data to deal with. We'll skip this
        % part of the data anyway for analysis, so no hurt incurred
        
        continue;
    end
    
    % value to hold during pursuit interval 
    if p==1 && pur.on(p)==1
        % use value after pursuit end as value to hold, since we have
        % nothign before
    else
        valX = X(pur.on(p)-1);
        valY = Y(pur.on(p)-1);
    end
    % do this carefully, see step from last pursuit sample to first saccade
    % sample and take this into account
    if pur.off(p) < nSamp
        stepX = diff(X(pur.off(p)+[0 1]));
        stepY = diff(Y(pur.off(p)+[0 1]));
        % now shift whole trace behind
        X(pur.off(p)+1:end) = X(pur.off(p)+1:end) - X(pur.off(p)) + valX + stepX;
        Y(pur.off(p)+1:end) = Y(pur.off(p)+1:end) - Y(pur.off(p)) + valY + stepY;
    end
    % now hold val during interval
    X(pur.on(p):pur.off(p)) = valX;
    Y(pur.on(p):pur.off(p)) = valY;
end

% as most Ss show a preference to make saccades in a certain direction,
% integrated position trace shows a trend. Remove trend by
% least-squares procedure. Expect this trend to be approximately opposite
% to trend of desaccaded pursuit of course.
% See Collewijn & Tamminga - 1984 - Human smooth and saccadic eye movements
% during voluntary pursuit of different target motions on different
% backgrounds. J Physiol. 351, pp. 217-250.
s = [1:length(X)].';
t = regstats(X,s,'linear',{'beta'});
X = X-s*t.beta(2);
t = regstats(Y,s,'linear',{'beta'});
Y = Y-s*t.beta(2);