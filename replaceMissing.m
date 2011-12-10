function data = replaceMissing(data,qInterpMissingPos)
% Currently, this only replaces missing data with linear interpolation for
% the velocity traces. It should be straightforward to implement if you
% need this for position or other traces.

% process data
[data.deg.Azi,...
 data.deg.Ele,...
 data.deg.vel,...
 data.deg.velAzi,...
 data.deg.velEle] = replaceMissingImplementation(data.deg.Azi,...
                                                 data.deg.Ele,...
                                                 data.deg.vel,...
                                                 data.deg.velAzi,...
                                                 data.deg.velEle,...
                                                 qInterpMissingPos,true);

if isfield(data.pix,'vel')
    [data.pix.X,...
     data.pix.Y,...
     data.pix.vel,...
     data.pix.velX,...
     data.pix.velY] = replaceMissingImplementation(data.pix.X,...
                                                   data.pix.Y,...
                                                   data.pix.vel,...
                                                   data.pix.velX,...
                                                   data.pix.velY,...
                                                   qInterpMissingPos,false);
end




function [X,Y,vel,velX,velY] = replaceMissingImplementation(X,Y,vel,velX,velY,qInterpMissingPos,qPrintInfo)

% We want to deal with all the nan in the data.
% This is getting rid of blinks and such...
qNaN = isnan(vel);
if any(qNaN)
    if qPrintInfo
        fprintf('  N NaN samples: %d (%.2f%%)\n',sum(qNaN),sum(qNaN)./length(vel)*100);
    end
    
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
    
    % adjust indices, nanon(p) and nanoff(p) point to first and last NaN in
    % a run
    nanon  = nanon-1;
    nanoff = nanoff+1;
    nanlen = nanoff-nanon+1;
    
    for p=1:length(nanon)
        % replace with interpolated velocity
        [vel,velX,velY] = replaceIntervalVelocity(vel,velX,velY,nanon(p),nanoff(p));
        
        if qInterpMissingPos
            X(nanon(p):nanoff(p)) = linspace(X(nanon(p)), X(nanoff(p)), nanlen(p));
            Y(nanon(p):nanoff(p)) = linspace(Y(nanon(p)), Y(nanoff(p)), nanlen(p));
        end
    end
    
    if qPrintInfo
        % show how many NaN we have left now, those cannot be handled
        fprintf('   -> N NaN samples left: %d (%.2f%%)\n',sum(isnan(vel)),sum(isnan(vel))./length(vel)*100);
    end
end