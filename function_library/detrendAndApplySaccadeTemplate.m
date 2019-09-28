function data = detrendAndApplySaccadeTemplate(data,ETparams)

% this function detrends the data (removes smooth pursuit) by first
% calculating the median velocity in a moving window, after some low-pass
% filtering to get an estimate of speed relatively unaffected by noise.
% This median velocity is then subtracted from the eye data.

% then, optionally, apply an FIR filter to this data to bring out the
% saccades from the noise. This filter is shaped after the velocity profile
% of a saccade and thus brings out saccade-like features in the trace by
% reducing the amplitude of other aspects (noise)!

% detrend/remove pursuit by subtracting median in a moving window
if ETparams.data.detrendWithMedianFilter
    %%% prepare algorithm parameters
    medianSamples   = ceil(ETparams.data.medianWindowLength/1000 * ETparams.samplingFreq);
    
    % decide which traces to detrend
    % we'll always want to do eye velocity in degrees as thats what the
    % saccade classification code runs on
    todo            = {'deg','vel','velDetrend'};
    
    if ETparams.data.detrendAll
        if ETparams.data.alsoStoreComponentDerivs
            % also azimuthal and elevational velocities
            todo        = [todo;...
                            {'deg','velAzi','velAziDetrend';...
                             'deg','velEle','velEleDetrend'}];
        end
        if ETparams.data.alsoStoreandDiffPixels
            % also pixels
            todo        = [todo;...
                            {'pix','vel','velDetrend'}];
            if ETparams.data.alsoStoreComponentDerivs
                % also X and Y speeds
                todo    = [todo;...
                            {'pix','velX','velXDetrend';...
                             'pix','velY','velYDetrend'}];
            end
        end
    end
        
    % do detrending -- don't have to deal with nans here unless the median
    % filter width is smaller than half the lowpass filter's number of taps
    for p=1:size(todo,1)
        median_velocities = medianFilter(...
            applyFilter(data.(todo{p,1}).(todo{p,2}),ETparams.data.detrendLowpassFIRCoeffs), ...    % lowpass filter first
            medianSamples);                                                                         % then median filter
        
        data.(todo{p,1}).(todo{p,3}) = data.(todo{p,1}).(todo{p,2}) - median_velocities;
    end
end

if ETparams.data.applySaccadeTemplate
    % convolute (do "pattern matching")
    
    % choose data
    if ETparams.data.detrendWithMedianFilter
        field = 'velDetrend';
    else
        field = 'vel';
    end
    thisdata = data.deg.(field);
    
    % replace nan with 0 so our filter responses don't get cut short by NaN
    % -> maximize the extent of the data in which we can find velocity
    % peaks. Important as high speeds are common next to missing data
    qNan                    = isnan(thisdata);
    thisdata(qNan)          = 0;
    
    % cross correlate, dealing with end effects by putting repeating first/last valid output
    correlation_responses   = applyFilter(thisdata,ETparams.data.saccadeTemplate);
    % scale and take absolute
    data.deg.velXCorr       = abs(correlation_responses .* ETparams.data.saccadeTemplateFilterScale);
                            
    % put nans back in as long runs of zero might bias the threshold
    % estimation step to lower thresholds.
    data.deg.velXCorr(qNan) = nan;
end