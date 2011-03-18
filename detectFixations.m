function data = detectFixations(data,ETparams)
%--------------------------------------------------------------------------
% Fixation detection
% Fixation are detected implicitly, as sections that are not saccade or
% glissade and also fit some other criteria
%--------------------------------------------------------------------------

%%% find data that is not saccade or glissade
[fixon,fixoff]  = bool2bounds(~(...
    bounds2bool(data.saccade.on ,data.saccade.off ,length(data.deg.vel)) | ...
    bounds2bool(data.glissade.on,data.glissade.off,length(data.deg.vel))   ...
    ));

%%% prepare algorithm parameters
minFixSamples   = ceil(ETparams.fixation.minDur/1000 * ETparams.samplingFreq);

%%% Process the tentative fixations
for kk = length(fixon):-1:1
    % Check that the fixation duration exceeds the minimum duration
    % criterion. Delete if not
    if fixoff(kk)-fixon(kk)+1 < minFixSamples
        fixon (kk) = [];
        fixoff(kk) = [];
        continue;
    end
    
    % Exclude section if any of the samples has a velocity > peak saccade
    % threshold, it cannot be a fixation (section somehow got deleted by
    % the saccade algorithm)
    if any(data.deg.vel(fixon(kk):fixoff(kk)) > data.saccade.peakVelocityThreshold)
        fixon (kk) = [];
        fixoff(kk) = [];
        continue;
    end
    
    % If the fixation contains NaN samples, handle it
    if any(isnan(data.deg.vel(fixon(kk):fixoff(kk))))
        switch ETparams.fixation.treatNaN
            
            case 1
                % NaN during fixation not allowed, so delete fixation
                fixon (kk) = [];
                fixoff(kk) = [];
                continue;
                
            case 2
                % allow NaN, simply ignore when computing average
                % fixation position (uses nanmean)
                % do nothing
                
            case 3
                % split fixation into multiple pieces, each without NaNs.
                % Keep pieces if they meet the minimum duration criterion
                % (we don't have to worry about the peak velocity
                % criterion). NB: The data markers in this piece of code
                % are local to the piece of data that is the current
                % "fixation" begin processed.
                
                % find non-NaN sections of data
                [dataon,dataoff] = bool2bounds(~isnan(data.deg.vel(fixon(kk):fixoff(kk))));
                
                % delete markers for the current fixation. we might
                % generate some new ones to replace it below
                % remember beginmark
                beginfixmark = fixon(kk);
                fixon (kk) = [];
                fixoff(kk) = [];
                
                % calculate length of each non-NaN section
                datalens = dataoff-dataon+1;
                
                % see if any pieces of data are long enough
                qLongEnough = datalens >= minFixSamples;
                
                if ~any(qLongEnough)
                    % if not, we're done
                    continue;
                end
                
                % if we have found some, isolate these good sections and
                % insert them in the fixation marker vectors
                dataon (~qLongEnough) = [];
                dataoff(~qLongEnough) = [];
                
                % markers are relative to start of this interval, add it
                dataon  = dataon  + beginfixmark-1;
                dataoff = dataoff + beginfixmark-1;
                
                % insert back
                fixon  = [ fixon(1:kk-1) dataon   fixon(kk:end)];
                fixoff = [fixoff(1:kk-1) dataoff fixoff(kk:end)];
        end
    end
    
    %%%%
    % Done. All the above criteria are fulfilled, we've got a fixation.
    %%%%
end

%%% output
data.fixation.on    = fixon;
data.fixation.off   = fixoff;