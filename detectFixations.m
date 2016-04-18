function data = detectFixations(data,ETparams)
%--------------------------------------------------------------------------
% Fixation detection
% Fixation are detected implicitly, as sections that are not saccade or
% glissade and also fit some other criteria
%--------------------------------------------------------------------------

%%% find data that is not saccade, glissade or blink
qNonFixFlag = bounds2bool(data.saccade .on,data.saccade .off, length(data.deg.vel));
if isfield(data,'glissade')
    qNonFixFlag = qNonFixFlag | bounds2bool(data.glissade.on,data.glissade.off, length(data.deg.vel));
end
if isfield(data,'blink')
    qNonFixFlag = qNonFixFlag | bounds2bool(data.blink   .on,data.blink   .off, length(data.deg.vel));
end
[fixon,fixoff]  = bool2bounds(~qNonFixFlag);

% correct so that fixation ends overlap with saccade starts (and etc)
% instead of 1 sample offset
fixon   = max(fixon -1, 1);
fixoff  = min(fixoff+1, length(data.deg.vel));

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
    % the saccade algorithm). Ignore if this happens at very end of data,
    % those velocity estimates are not so reliable
    if any(data.deg.vel(max(3,fixon(kk)):min(fixoff(kk),length(data.deg.vel)-2)) > data.saccade.peakVelocityThreshold)
        fixon (kk) = [];
        fixoff(kk) = [];
        continue;
    end
    
    % If the fixation contains NaN samples, handle it
    qNaN = isnan(data.deg.vel(fixon(kk):fixoff(kk)));
    if any(qNaN)
        if all(qNaN)
            % no fixation data, remove
            fixon (kk) = [];
            fixoff(kk) = [];
            continue;
        end
        switch ETparams.fixation.treatNaN
            
            case 1
                % NaN during fixation not allowed, so delete fixation
                fixon (kk) = [];
                fixoff(kk) = [];
                continue;
                
            case {2,3}
                % 2: allow NaN, simply ignore when computing average
                % fixation position (uses nanmean). However, make sure eye
                % position doesn't jump significantly during the missing
                % data.
                % 3: split fixation into multiple pieces, each without
                % NaNs. Keep pieces if they meet the minimum duration
                % criterion (we don't have to worry about the peak velocity
                % criterion).
                % NB: The data markers in this piece of code are local to
                % the piece of data that is the current "fixation" begin
                % processed.
                
                % find non-NaN sections of data
                [dataon,dataoff] = bool2bounds(~qNaN);
                
                % delete markers for the current fixation. we might
                % generate some new ones to replace it below
                % remember beginmark
                beginfixmark = fixon(kk);
                fixon (kk) = [];
                fixoff(kk) = [];
                
                if ETparams.fixation.treatNaN==2
                    startidxs   = beginfixmark+dataoff(1:end-1)-1;
                    endidxs     = beginfixmark+dataon (2:end)  -1;
                    
                    % calculate eye position change between start and end
                    % of nan sections.
                    amps        = calcAmplitudeFick(...
                                      data.deg.Azi(startidxs), data.deg.Ele(startidxs),...
                                      data.deg.Azi(endidxs  ), data.deg.Ele(endidxs  )...
                                  );
                    
                    okIdx       = find(amps<ETparams.fixation.NaNMaxJump);
                    
                    % remove jumps that are not over maximum amplitude from
                    % processing queue of possible split points
                    dataon (okIdx+1) = [];
                    dataoff(okIdx  ) = [];
                end
                
                % calculate length of each (non-NaN) section
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
