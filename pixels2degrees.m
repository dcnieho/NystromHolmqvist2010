function [alphaH,alphaV] = pixels2degrees(pixels)
% Converts pixels to degrees
% Input: 
global ETparams

d = ETparams.viewingDist;
pixelDimensions = ETparams.screenSz;
stimulusDimensions = ETparams.screenDim;

meterPerPixel = pixels*stimulusDimensions./pixelDimensions;

alphaH = 180/pi*(2*atan(meterPerPixel(1)/d/2));
alphaV = 180/pi*(2*atan(meterPerPixel(2)/d/2));

