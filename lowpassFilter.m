function [output] = lowpassFilter(input,filter_coeffs)

% function [output] = lowpass_filter(input);
%


if 0
    % old code:
    % apply the 2-d FIR filter described by filter_coeffs
    output2 = filter2(filter_coeffs,input);
    
    % chad padded one filter lengths' worth of data
    % at the outset and end of each data channel
    output2(1:length(filter_coeffs)) = output2(length(filter_coeffs)+1);
    output2(end-length(filter_coeffs):end) = output2(end-length(filter_coeffs));
else
    % NB, note the same as the above! The above dealt with edge effects too
    % rigorously, we don't have to add a whole filter's length, just
    % half
    ntap    = length(filter_coeffs);
    window  = floor(ntap/2);
    % apply the 2-d FIR filter described by filter_coeffs
    output  = filter2(filter_coeffs,input,'valid');
    % deal with end effects
    output  = [repmat(output(1),window,1); output; repmat(output(end),window,1)];
end