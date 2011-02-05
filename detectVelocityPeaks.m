function detectVelocityPeaks(i,j)
global ETparams
 

 % Find peaks larger than a threshold ('ETparams.peakDetectionThreshold')
 % Sets a '1' where the velocity is larger than the threshold and '0'
 % otherwise
ETparams.data(i,j).InitialVelPeakIdx  = (ETparams.data(i,j).vel > ETparams.data(i,j).peakDetectionThreshold);
 
 
