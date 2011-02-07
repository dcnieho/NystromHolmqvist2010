function eventDetection
global ETparams

%--------------------------------------------------------------------------
% Process eye-movement data for file (participant) and trial separately, i - files, j -
% trials
%--------------------------------------------------------------------------
fprintf('%s\n','Detecting events')
for i = 2%:size(ETparams.data,1)

    for j = 1:size(ETparams.data,2)
        data = ETparams.data(i,j);
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
        % temp: make figure, TODO: move this to separate functions soon
        figure,plot(data.vel,'k')
        hold on
        plot([0 length(data.vel)],[1 1]*data.peakDetectionThreshold,'r--')
        plot([0 length(data.vel)],[1 1]*data.saccadeVelocityTreshold,'r:')
        plot(data.saccade.on,data.vel(data.saccade.on),'bo','MarkerFaceColor','blue','MarkerSize',4)
        plot(data.saccade.off,data.vel(data.saccade.off),'ro','MarkerFaceColor','red','MarkerSize',4)
        qhighvelglissade = data.glissade.type==2;
        plot(data.glissade.off(~qhighvelglissade),data.vel(data.glissade.off(~qhighvelglissade)),'g*')
        plot(data.glissade.off( qhighvelglissade),data.vel(data.glissade.off( qhighvelglissade)),'c*')
        for p=1:length(data.saccade.off)
            plot([0 ceil(ETparams.glissadeSearchWindowms/1000 * ETparams.samplingFreq)]+data.saccade.off(p),[1 1]*data.saccade.localVelocityTreshold(p),'r-');
        end
        % end temp

        % Implicitly detect fixations
        %-------------------------------------            
        data = detectFixations(data,ETparams);
    end
end