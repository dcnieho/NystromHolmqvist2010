function [median_channel] = medianFilter(input,window_length)

array           = rovingBin(input,window_length);
median_channel  = median(array,2,'omitnan');