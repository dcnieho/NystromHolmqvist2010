function data = processSaccadesAndGlissades(data,ETparams)
% gets information about saccades and glissades

% Collect information about the saccade
data    = processSaccadesAndGlissadesImplementation(data, 'saccade' ,ETparams);

% Collect information about the glissade
if isfield(data,'glissade')
    data    = processSaccadesAndGlissadesImplementation(data, 'glissade',ETparams);
end






% helper function (as we need to run exactly the same code over saccades
% and glissades)
function data = processSaccadesAndGlissadesImplementation(data, fieldname, ETparams)

% duration (in milliseconds)
data.(fieldname).duration = (data.(fieldname).off-data.(fieldname).on+1)/ETparams.samplingFreq;

% amplitude & direction
[data.(fieldname).amplitude, data.(fieldname).direction] = ...
    calcAmplitudeFick(...
        data.deg.Azi(data.(fieldname).on ), data.deg.Ele(data.(fieldname).on ),...
        data.deg.Azi(data.(fieldname).off), data.deg.Ele(data.(fieldname).off)...
    );

% and now some that are best done in a for-loop
data.(fieldname).peakVelocity     = zeros(size(data.(fieldname).on));
data.(fieldname).peakAcceleration = zeros(size(data.(fieldname).on));
for p=1:length(data.(fieldname).on)
    idxs = data.(fieldname).on(p) : data.(fieldname).off(p);
    data.(fieldname).peakVelocity(p)       = max(data.deg.vel(idxs));
    data.(fieldname).peakAcceleration(p)   = max(data.deg.acc(idxs));
end