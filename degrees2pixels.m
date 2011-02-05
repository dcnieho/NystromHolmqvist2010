function [angleInPixelsH, angleInPixelsV] = degrees2pixels(alpha,d,pixelDimensions,stimulusDimensions)

alpha = alpha*pi/180;

width05 = tan(alpha/2)*d;
width = 2*width05;

angleInPixelsH = pixelDimensions(1)*width/stimulusDimensions(1);
angleInPixelsV = pixelDimensions(2)*width/stimulusDimensions(2);