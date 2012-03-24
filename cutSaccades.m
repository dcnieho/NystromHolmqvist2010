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
% deal with velocity
for p=1:length(sac.on) * ~~bitand(cutPostrace,uint8(1))
    % per saccade, linearly interpolate eye velocity between saccade begin
    % point and end point
    if ~~bitand(cutPostrace,uint8(8)) && sac.on(p) <= skipWindowSamples
        % special case: when reconstructing eye position, skip saccades
        % during the first second, they are probably to catch the blob as
        % it appears at the beginning of the trial. Otherwise we'll have a
        % constant offest in the data to deal with. We'll skip this part of
        % the data anyway for analysis, so no hurt incurred
        
        continue;
    end
    
    % replace with interpolated velocity
    [vel,velX,velY] = replaceIntervalVelocity(vel,velX,velY,sac.on(p),sac.off(p));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with position
switch bitand(cutPostrace,uint8(4+8+16))
    case 8
        % integrate velocities using plant
        X = CanonicalDiscreteSSModel(ETparams.sysdt,velX).' + X(1);
        Y = CanonicalDiscreteSSModel(ETparams.sysdt,velY).' + Y(1);
    case 16
        % simply linearly interpolate position
        for p=1:length(sac.on)
            X(sac.on(p):sac.off(p)) = linspace(X(sac.on(p)), X(sac.off(p)), sac.len(p));
            Y(sac.on(p):sac.off(p)) = linspace(Y(sac.on(p)), Y(sac.off(p)), sac.len(p));
        end
    case 4
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
                % zorg at we niet buiten de data grijpen
                start1 = max(start1,1);
            else
                % zorg at we niet een stuk van een vorige saccade meenemen
                start1 = max(start1,sac.off(p-1));
            end
            eind1  = sac.on(p)-eoff;
            
            start2 = sac.off(p)+eoff;
            eind2  = sac.off(p)+soff;
            if p==length(sac.on)
                % zorg at we niet buiten de data grijpen
                eind2 = min(eind2,length(X));
            else
                % zorg at we niet een stuk van een volgende saccade meenemen
                eind2 = max(eind2,sac.on(p+1));
            end
            
            X = interpSacIgnace(X,2,start1,eind1,start2,eind2,[sac.on(p) sac.off(p)]); % add ,true); fr debug plots
            Y = interpSacIgnace(Y,2,start1,eind1,start2,eind2,[sac.on(p) sac.off(p)]);
        
            
            % now, reconstruct velocity from desaccaded position, if wanted.
            % calc velocity during interpolated saccades and put that in
            if ~~bitand(cutPostrace,uint8(2))
                on  = sac.on (p);
                off = sac.off(p);
                
                velX(on:off) = diff(X(on:off+1)) * ETparams.samplingFreq;
                velY(on:off) = diff(Y(on:off+1)) * ETparams.samplingFreq;
                
                vel(on:off)  = hypot(velX(on:off),velY(on:off));
            end
        end
end