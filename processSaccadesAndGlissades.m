function data = processSaccadesAndGlissades(data,ETparams)
% merges and then gets information about saccades and glissades

%--------------------------------------------------------------------------
% MERGING SACCADES
%--------------------------------------------------------------------------
% now deal with saccades that are too close together, fuse two saccades
% with little time between them. This needs to run even if the window
% length is 0 as the above has been seen to generate saccade starts before
% the end of the previous saccade. We could check for that above, or we
% could just prune/merge them here.

% prepare parameters
% if saccade is followed by another saccade by less than SacMergeWindow (in
% ms), they'll be merged
SacMergeWindowSamp = ceil(ETparams.saccade.mergeWindow./1000 * ETparams.samplingFreq);
if ETparams.data.qApplySaccadeTemplate && ETparams.saccade.qSaccadeTemplateRefine
    field_offset= 'offsetXCorrThreshold';
else
    field_offset= 'offsetVelocityThreshold';
end

kk=1;
while kk < length(data.saccade.on)      % NB: doesn't process last saccade (useless anyway of course!)
    % walk through all saccades and see if followed shortly by another
    % saccade.
    
    % If there is a glissades for this saccade, check from the end of the
    % glissade, otherwise check from the end of the saccade
    if ETparams.glissade.qDetect
        qHaveGlissade = data.saccade.off(kk)==data.glissade.on;
    else
        qHaveGlissade = false;
    end
    if any(qHaveGlissade)
        assert(sum(qHaveGlissade)==1)  % anything else would be ridiculous!
        thisoff = data.glissade.off(qHaveGlissade);
    else
        thisoff = data.saccade.off(kk);
    end
    
    % check if start of next saccade occurs within the window
    if data.saccade.on(kk+1)-thisoff <= SacMergeWindowSamp
        % if yes, merge
        data.saccade.off(kk)                = data.saccade.off(kk+1);
        
        % remove next saccade...
        data.saccade.on(kk+1)               = [];
        data.saccade.off(kk+1)              = [];
        data.saccade.(field_offset)(kk+1)   = [];
        
        % ... and glissade that is caught in between
        if qHaveGlissade
            data.glissade.on (qHaveGlissade)    = [];
            data.glissade.off(qHaveGlissade)    = [];
        end
        
        continue;
    else
        % else, process next saccade
        kk = kk+1;
    end
end


%--------------------------------------------------------------------------
% COLLECT INFORMATION
%--------------------------------------------------------------------------

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