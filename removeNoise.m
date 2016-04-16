function data = removeNoise(data,ETparams)
% find those sections of data enclosed in nans that are too short to be
% meaningful (less than minimum duration provided by user) and delete them

% prepare algorithm parameters
minSamples   = ceil(ETparams.data.minDur/1000 * ETparams.samplingFreq);

% process to look for short islands of data in the sea of missingness
[dataon,dataoff] = bool2bounds(~isnan(data.deg.vel));
for p=length(dataon):-1:1
    % Check that the section of data is longer than the minimum duration.
    % Keep the indices if not so we can delete it later
    if dataoff(p)-dataon(p)+1 >= minSamples
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
qNaN        = isnan(data.deg.vel);
data.deg    = replaceElementsInStruct(data.deg,qNaN,nan);

% if we have smoothed eye position in pixels and its derivatives, throw the
% NaNs in there as well
if isfield(data.pix,'X')
    data.pix    = replaceElementsInStruct(data.pix,qNaN,nan);
end

% lastly, notify if more than 20% nan
if sum(qNaN)/length(data.deg.vel) > 0.20
    fprintf('Warning: This trial contains %.2f%% missing samples\n',sum(qNaN)/length(data.deg.vel)*100);
    data.qNoiseTrial = true;
else
    data.qNoiseTrial = false;
end