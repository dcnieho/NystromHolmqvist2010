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
if ETparams.data.qDetrendWithMedianFilter
    %%% prepare algorithm parameters
    medianSamples   = ceil(ETparams.data.medianWindowLength/1000 * ETparams.samplingFreq);
    
    % decide which traces to detrend
    % we'll always want to do eye velocity in degrees as thats what the
    % saccade detection code runs on
    todo            = {'deg','vel','velDetrend'};
    
    if ETparams.data.qDetrendAll
        if ETparams.data.qAlsoStoreComponentDerivs
            % also azimuthal and elevational velocities
            todo        = [todo;...
                            {'deg','velAzi','velAziDetrend';...
                             'deg','velEle','velEleDetrend'}];
        end
        if ETparams.data.qAlsoStoreandDiffPixels
            % also pixels
            todo        = [todo;...
                            {'pix','vel','velDetrend'}];
            if ETparams.data.qAlsoStoreComponentDerivs
                % also X and Y speeds
                todo    = [todo;...
                            {'pix','velX','velXDetrend';...
                             'pix','velY','velYDetrend'}];
            end
        end
    end
        
    % do detrending
    for p=1:size(todo,1)
        median_velocities = medianFilter(...
            lowpassFilter(data.(todo{p,1}).(todo{p,2}),ETparams.data.detrendLowpassFIRCoeffs), ...
            medianSamples); % smoothed median matches
        
        data.(todo{p,1}).(todo{p,3}) = data.(todo{p,1}).(todo{p,2}) - median_velocities;
    end
end

if ETparams.data.qApplySaccadeTemplate
    % convolute and deal with the end effects
    % part of convolution where filter runs off the data, at both ends
    trans                   = floor(length(ETparams.data.saccadeTemplate)/2);
    % choose data
    if ETparams.data.qDetrendWithMedianFilter
        field = 'velDetrend';
    else
        field = 'vel';
    end
    thisdata = data.deg.(field);
    
    % replace nan with 0 so our filter responses don't get cut short by NaN
    % -> maximize the extent of the data in which we can detect velocity
    % peaks. Important as high speeds are common next to missing data
    qNan = isnan(thisdata);
    thisdata(qNan) = 0;
    
    % cross correlate, dealing with end effects by putting NaN there
    correlation_responses   = [NaN(trans,1); ...
                                conv(thisdata,ETparams.data.saccadeTemplate,'valid'); ...
                                NaN(trans,1)];
                            
    % put nans back in as long runs of zero might bias the threshold
    % estimation step to lower thresholds. Only remove where filter would
    % have been running on nan only
    [nanon,nanoff]  = bool2bounds(qNan);
    nanon           = nanon +trans;
    nanoff          = nanoff-trans;
    correlation_responses(bounds2ind(nanon,nanoff)) = nan;
    
    % scale and take absolute
    data.deg.velXCorr       = abs(correlation_responses .* ETparams.data.saccadeTemplateFilterScale);
end