function ETparams = prepareParameters(ETparams)

% calculate extends of screen around the origin (we'll change the data's
% origin in the same way in prepareData)
ETparams.screen.rect.pix = [...
     [0 0]                      - ETparams.screen.subjectStraightAhead, ...
     ETparams.screen.resolution - ETparams.screen.subjectStraightAhead
    ];

% calculate screen extends in degree
pixPerMeter     = ETparams.screen.resolution ./ ETparams.screen.size;
screenExtends   = ETparams.screen.rect.pix ./ [pixPerMeter pixPerMeter];
ETparams.screen.rect.deg = atan2(screenExtends,ETparams.screen.viewingDist)*180/pi;

% calculate low-pass filter for median filter detrending step, if needed
if ETparams.data.qDetrendWithMedianFilter
    % -3 dB rolloff at 20Hz, using 7 coefficient convolution
    % (at 240Hz)
    %
    % filter_coeffs_240 = [.0303 .1639 .2359 .2838 .2359 .1639 .0303]
    
    filter_coeffs_240 = [.0303 .1639 .2359 .2838 .2359 .1639 .0303];
    
    if ETparams.samplingFreq ~= 240
        % resample at ETparams.samplingFreq
        p = polyfit(-3:3,filter_coeffs_240,4);  % fourth order polynomial seems to interpolate this one very nicely!
        num_samples = round(ETparams.samplingFreq/240 * length(filter_coeffs_240));
        filter_coeffs = polyval(p,linspace(-3,3,num_samples));
    else
        filter_coeffs = filter_coeffs_240;
    end
    
    % normalize
    filter_coeffs = filter_coeffs ./ sum(filter_coeffs);
    
    % ensure column vector
    ETparams.data.detrendLowpassFIRCoeffs = filter_coeffs(:);
    
    if 0
        % Debug: plot resampled filter coefficients
        figure(200)
        clf
        plot(-3:3,filter_coeffs_240,'.');
        hold on
        
        plot(-3:.01:3,polyval(p,-3:.01:3),'g')
        plot(linspace(-3,3,num_samples),ETparams.data.detrendLowpassFIRCoeffs,'r.');
    end
end

% calculate saccade template FIR
if ETparams.data.qApplySaccadeTemplate
    lstone_velocity_saccade_template_7 = [.03 .106 .221 .285 .221 .106 .03];
    saccade_template_240 = lstone_velocity_saccade_template_7;
    
    if ETparams.samplingFreq ~= 240
        % resample at sampling_rate
        p = polyfit(-3:3,saccade_template_240,6);  % sixth order polynomial seems better for this one
        num_samples = round(ETparams.samplingFreq/240 * length(saccade_template_240));
        saccade_template = polyval(p,linspace(-3,3,num_samples));
    else
        saccade_template = saccade_template_240;
    end
    
    % normalize template
    if 1
        saccade_template = saccade_template ./ sum(saccade_template);
    end
    
    % ensure column vector
    ETparams.data.saccadeTemplate = saccade_template(:);
    
    if 0
        % plot resampled saccade template
        figure(201)
        clf
        plot(-3:3,saccade_template_240,'.');
        hold on
        
        plot(-3:.01:3,polyval(p,-3:.01:3),'g')
        plot(linspace(-3,3,num_samples),ETparams.data.saccadeTemplate,'r.');
    end
    
    % calculate filter scale for normalization
    saccade_size = sum(ETparams.data.saccadeTemplate) / ETparams.samplingFreq;
    template_SS  = sum(ETparams.data.saccadeTemplate.^2);
    ETparams.data.saccadeTemplateFilterScale = saccade_size ./ template_SS;
end