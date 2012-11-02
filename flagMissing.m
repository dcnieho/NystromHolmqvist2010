function data = flagMissing(data)

[data.deg.missing] = flagMissingImplementation(data.deg.vel,true);

if isfield(data.pix,'vel')
    [data.pix.missing] = flagMissingImplementation(data.pix.vel,false);
end




function [missing] = flagMissingImplementation(vel,qPrintInfo)


missing.on  = [];
missing.off = [];

% We want to deal with all the nan in the data.
% This is getting rid of blinks and such...
qNaN = isnan(vel);
if any(qNaN)
    if qPrintInfo
        fprintf('  N NaN samples: %d (%.2f%%)\n',sum(qNaN),sum(qNaN)./length(vel)*100);
    end
    
    [nanon,nanoff] = bool2bounds(qNaN);
    missing.on  = nanon;
    missing.off = nanoff;
end