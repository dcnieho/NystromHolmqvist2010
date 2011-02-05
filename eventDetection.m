function eventDetection
global ETparams

%--------------------------------------------------------------------------
% Process eye-movement data for file (participant) and trial separately, i - files, j -
% trials
%--------------------------------------------------------------------------
fprintf('%s\n','Detecting events')
for i = 2%:size(ETparams.data,1)

    for j = 1:size(ETparams.data,2)
%         i,j
        % Calculate velocity and acceleration
        %-------------------------------------
        data = calVelAcc_sgolay(ETparams.data(i,j),ETparams);
        % temporary for compatibility
        fields = fieldnames(data);
        for p=1:length(fields)
            ETparams.data(i,j).(fields{p}) = data.(fields{p});
        end
        % end temp

        % Detect blinks and noise
        %-------------------------------------
        [data,qnoise] = detectAndRemoveNoise(ETparams.data(i,j),ETparams);
        % temporary for compatibility
        fields = fieldnames(data);
        for p=1:length(fields)
            ETparams.data(i,j).(fields{p}) = data.(fields{p});
        end
        ETparams.nanIdx(i,j).Idx = qnoise;
        % end temp

        % iteratively find the optimal noise threshold
        %-------------------------------------
        
        data = estimateSaccadeVelocityThresholds(ETparams.data(i,j),ETparams);
        % temporary for compatibility
        fields = fieldnames(data);
        for p=1:length(fields)
            ETparams.data(i,j).(fields{p}) = data.(fields{p});
        end
        % end temp
        
        % Detect saccades (with peak detection threshold (v < v_avg_noise + 3*v_std_noise))
        % and glissades
        %-------------------------------------            
        detectSaccades(i,j)

        % Implicitly detect fixations
        %-------------------------------------            
        detectFixations(i,j)

    end
end

    

