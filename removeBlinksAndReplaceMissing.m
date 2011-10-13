function data = removeBlinksAndReplaceMissing(data)
% Currently, this only replaces missing data with linear interpolation for
% the velocity traces. It should be straightforward to implement if you
% need this for position or other traces.
% Missing data during saccades is flagged as blinks. With the Eyelink,
% blinks lead to rapid downward/upward eyemovements being detected
% preceding and following it. These saccades are thrown out.
% NB: You would need to rerun processSaccadesAndGlissades.m after this as
% the saccades flags have just changed

[data,qBlinkDeg] = flagBlinksAndRemoveMissing(data,'deg',true);

if isfield(data.pix,'vel')
    [data,qBlinkPix] = flagBlinksAndRemoveMissing(data,'pix',false);
    
    % there should be no difference in missing between the pixel and degree
    % traces so this should always be true    assert(isequal(qBlinkPix,qBlinkDeg))
end

% now remove saccades that were flagged as blink
sacFields   = fieldnames(data.saccade);
qScalar     = structfun(@isscalar,data.saccade);
sacFields(qScalar) = [];
for p=1:length(sacFields)
    data.saccade.(sacFields{p})(qBlinkDeg) = [];
end



function [data,qBlink] = flagBlinksAndRemoveMissing(data,datatype,qPrintInfo)

% get eye velocities in pixels/degree
vel     = data.(datatype).vel;
if strcmp(datatype,'pix')
    velX    = data.pix.velX;
    velY    = data.pix.velY;
elseif strcmp(datatype,'deg')
    velX    = data.deg.velAzi;
    velY    = data.deg.velEle;
end

sacon       = data.saccade.on;
sacoff      = data.saccade.off;

qBlink      = false(size(sacon));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step one: we want to deal with all the nan in the data.
% This is getting rid of blinks and such...
qNaN = isnan(vel);
if any(qNaN)
    if qPrintInfo
        fprintf('  N NaN samples: %d\n',sum(qNaN));
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
    
    for p=1:length(nanon)
        qDuringSac = nanon(p)>=sacon & nanoff(p)<=sacoff;
        if any(qDuringSac)
            % if nan is during saccade, use those as start and end points
            % (probably blink)
            assert(sum(qDuringSac)==1)  % anything else would be ridiculous!
            on  = sacon (qDuringSac);
            off = sacoff(qDuringSac);
            
            % flag blink for later removal
            qBlink = qBlink | qDuringSac;
        else
            % pas indices aan, nanon(p) and nanoff(p) wijzen naar de eerste
            % en laatste NaN in een serie
            on  = nanon(p)-1;
            off = nanoff(p)+1;
        end
        % replace with interpolated velocity
        [vel,velX,velY] = replaceIntervalVelocity(vel,velX,velY,on,off);
    end
    
    if qPrintInfo
        % show how many NaN we have left now, those cannot be handled
        fprintf('   -> N NaN samples left: %d\n',sum(isnan(vel)));
    end
    
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