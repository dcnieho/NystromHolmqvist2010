function data = eventDetection(data,ETparams)

% Calculate velocity and acceleration
%-------------------------------------
data = calVelAcc_sgolay(data,ETparams);

% Detect blinks and noise
%-------------------------------------
data= detectAndRemoveNoise(data,ETparams);

% iteratively find the optimal noise threshold
%-------------------------------------
data = estimateSaccadeVelocityThresholds(data,ETparams);

% Detect saccades (with peak detection threshold (v < v_avg_noise + 3*v_std_noise))
% and glissades
%-------------------------------------
data = detectSaccades(data,ETparams);

% Implicitly detect fixations
%-------------------------------------
data = detectFixations(data,ETparams);