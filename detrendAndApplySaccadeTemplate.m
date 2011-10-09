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
    todo            = {'deg','vel','detrend_vel'};
    
    if ETparams.data.qDetrendAll
        if ETparams.data.qAlsoStoreComponentDerivs
            % also azimuthal and elevational velocities
            todo        = [todo;...
                            {'deg','velAz','detrend_velAz';...
                             'deg','velEl','detrend_velEl'}];
        end
        if ETparams.data.qAlsoStoreandSmoothPixels
            % also pixels
            todo        = [todo;...
                            {'pix','vel','detrend_vel'}];
            if ETparams.data.qAlsoStoreComponentDerivs
                % also X and Y speeds
                todo    = [todo;...
                            {'pix','velX','detrend_velX';...
                             'pix','velY','detrend_velY'}];
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
        field = 'detrend_vel';
    else
        field = 'vel';
    end
    
    % cross correlate, dealing with end effects by putting NaN there
    correlation_responses   = [NaN(trans,1); ...
                                conv(data.deg.(field),ETparams.data.saccadeTemplate,'valid'); ...
                                NaN(trans,1)];
    
    % scale and take absolute
    data.deg.xcorr_vel      = abs(correlation_responses .* ETparams.data.saccadeTemplateFilterScale);
end