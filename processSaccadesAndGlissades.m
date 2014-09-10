function data = processSaccadesAndGlissades(data,ETparams)
% merges and then gets information about saccades and glissades

%--------------------------------------------------------------------------
% MERGING SACCADES
%--------------------------------------------------------------------------
% now deal with saccades that are too close together, fuse two saccades
% with little time between them.

if ETparams.saccade.mergeWindow>0
    % prepare parameters
    % if saccade is followed by another saccade by less than mergeWindow (in
    % ms), they'll be merged
    SacMergeWindowSamples = ceil(ETparams.saccade.mergeWindow./1000 * ETparams.samplingFreq);
    
    if isfield(data,'glissade')
        [data.saccade,data.glissade] = mergeIntervals(data.saccade,data.glissade,SacMergeWindowSamples);
    else
        data.saccade                 = mergeIntervals(data.saccade,     []      ,SacMergeWindowSamples);
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

% start and end points
nSamp= ETparams.saccade.seWindowSamp;
sidx = bsxfun(@plus,[-nSamp:0    ].',data.(fieldname).on ); sidx(sidx<1) = 1;
data.(fieldname).startPointAzi = mean(data.deg.Azi(sidx),1);
data.(fieldname).startPointEle = mean(data.deg.Ele(sidx),1);
eidx = bsxfun(@plus,[     0:nSamp].',data.(fieldname).off); eidx(eidx>length(data.deg.Azi)) = length(data.deg.Azi);
data.(fieldname).endPointAzi   = mean(data.deg.Azi(eidx),1);
data.(fieldname).endPointEle   = mean(data.deg.Ele(eidx),1);

% use start and end points to calculate amplitude & direction
[data.(fieldname).amplitude, data.(fieldname).direction] = ...
    calcAmplitudeFick(...
        data.(fieldname).startPointAzi, data.(fieldname).startPointEle,...
        data.(fieldname).endPointAzi  , data.(fieldname).endPointEle...
    );

% and now some that are best done in a for-loop
data.(fieldname).peakVelocity     = zeros(size(data.(fieldname).on));
data.(fieldname).peakAcceleration = zeros(size(data.(fieldname).on));
for p=1:length(data.(fieldname).on)
    idxs = data.(fieldname).on(p) : data.(fieldname).off(p);
    data.(fieldname).peakVelocity(p)       = max(data.deg.vel(idxs));
    data.(fieldname).peakAcceleration(p)   = max(data.deg.acc(idxs));
end
