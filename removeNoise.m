function data = removeNoise(data,ETparams)
% find those sections of data enclosed in nans that are too short to be
% meaningful (less than minimum fixation duration) and delete them

% prepare algorithm parameters
minFixSamples   = ceil(ETparams.fixation.minDur/1000 * ETparams.samplingFreq);

% process
[dataon,dataoff] = bool2bounds(~isnan(data.deg.vel));
for p=length(dataon):-1:1
    % Check that the section of data is longer than the minimum fixation
    % duration. Keep the indices if not so we can delete it later
    if dataoff(p)-dataon(p) >= minFixSamples
        dataon (p) = [];
        dataoff(p) = [];
        continue;
    end
end
if ~isempty(dataon)
    noiseidxs = bounds2ind(dataon,dataoff);
    % remove useless data
    data.deg.vel(noiseidxs) = nan;
end

% now remove data from all other data fields that we have
qNaN = isnan(data.deg.vel);
fields = fieldnames(data.deg);
for p=1:length(fields)
    if strcmp(fields{p},'vel')
        continue;
    end
    
    data.deg.(fields{p})(qNaN) = nan;
end

% if we have smoothed eye position in pixels and its derivatives, throw the
% NaNs in there as well
if isfield(data.pix,'X')
    fields = fieldnames(data.pix);
    for p=1:length(fields)
        data.pix.(fields{p})(qNaN) = nan;
    end
end

% lastly, notify if more than 20% nan
if sum(isnan(data.deg.vel))/length(data.deg.vel) > 0.20
    fprintf('Warning: This trial contains %.2f%% missing samples\n',sum(isnan(data.deg.vel))/length(data.deg.vel)*100);
    data.qNoiseTrial = true;
else
    data.qNoiseTrial = false;
end