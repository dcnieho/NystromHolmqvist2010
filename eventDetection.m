function eventDetection
global ETparams

%--------------------------------------------------------------------------
% Process eye-movement data for file (participant) and trial separately, i - files, j -
% trials
%--------------------------------------------------------------------------
fprintf('%s','Detecting events')
for i = 1%:size(ETparams.data,1)

    for j = 1:size(ETparams.data,2)
%         i,j
        % Calculate velocity and acceleration
        %-------------------------------------
        if ETparams.qUseDN
            calVelAcc_sgolay_DN(i,j)
        else
            calVelAcc_sgolay(i,j)
        end

        % Detect blinks and noise
        %-------------------------------------
        detectAndRemoveNoise(i,j)

        % iteratively find the optimal noise threshold
        %-------------------------------------
        ETparams.data(i,j).peakDetectionThreshold = ETparams.peakDetectionThreshold;
        oldPeakT = inf;
        while abs(ETparams.data(i,j).peakDetectionThreshold -  oldPeakT) > 1

            oldPeakT  = ETparams.data(i,j).peakDetectionThreshold;

            % Detect peaks in velocity (> X degrees/second)
            detectVelocityPeaks(i,j) 

            % Find fixation noise level (0.7*global fixation noise +
            % 0.3*local fixation)
            detectFixationNoiseLevel(i,j)    

        end

        % Detect saccades (with peak detection threshold (v < v_avg_noise + 3*v_std_noise))
        % and glissades
        %-------------------------------------            
        detectSaccades(i,j)

        % Implicitly detect fixations
        %-------------------------------------            
        detectFixations(i,j)

    end
end

    

