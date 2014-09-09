function [median_channel] = medianFilter(input,window_length)

array   = rovingBin(input,window_length);

if 0
    % if no nanmedian from the statistics toolbox:
    median_channel = zeros(size(array,1),1);
    
    % KLUDGE that should be resolved with NANmedian, or similar
    for i=1:size(array,1)
        median_channel(i) = median(array(i,~isnan(array(i,:))));
    end
else
    median_channel = nanmedian(array,2);
end