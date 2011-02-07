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
        % temp: original code
        ETparams.data(i,j).velPeakIdx  = ETparams.data(i,j).vel > ETparams.data(i,j).peakDetectionThreshold;
        detectSaccadesOld(i,j);
        figure; hold on;
        plot(data.vel,'k')
        for kk = 1:size(ETparams.saccadeInfo,3)
            idxSaccadeEnd = ETparams.saccadeInfo(i,j,kk).end*ETparams.samplingFreq;
            idxSaccadeStart = ETparams.saccadeInfo(i,j,kk).start*ETparams.samplingFreq;
            
            if ~isempty(idxSaccadeEnd )
                plot(idxSaccadeStart,data.vel(round(idxSaccadeStart)),'bo');
                plot(idxSaccadeEnd,data.vel(round(idxSaccadeEnd)),'ro');
                % Plot local threshold
                plot((idxSaccadeStart:idxSaccadeEnd)*convf,ones(1,length(idxSaccadeStart:idxSaccadeEnd))*ETparams.data(i,j,kk).localSaccadeVelocityTreshold,'-k');
            end
        end
        
        % Plot glissades
        for kk = 1:size(ETparams.glissadeInfo,3)
            idxGlissadeEnd = ETparams.saccadeInfo(i,j,kk).end*ETparams.samplingFreq +  ETparams.glissadeInfo(i,j,kk).duration*ETparams.samplingFreq;
            if ETparams.glissadeInfo(i,j,kk).type==1
                plot(idxGlissadeEnd,data.vel(round(idxGlissadeEnd)),'g*');
            elseif ETparams.glissadeInfo(i,j,kk).type==2
                plot(idxGlissadeEnd,data.vel(round(idxGlissadeEnd)),'c*');
            end
        end
        for kk=1:size(ETparams.data,3)
            idxSaccadeEnd = ETparams.saccadeInfo(i,j,kk).end*ETparams.samplingFreq;
            if ~isempty(idxSaccadeEnd )
                plot([0 ceil(ETparams.glissadeSearchWindowms/1000 * ETparams.samplingFreq)]+idxSaccadeEnd,[1 1]*ETparams.data(i,j,p).localSaccadeVelocityTreshold,'--');
            end
        end
        ax = gca;
        % end temp
        data = detectSaccades(ETparams.data(i,j),ETparams);
        % temp: make figure
        figure,plot(data.vel,'k')
        hold on
        plot([0 length(data.vel)],[1 1]*data.peakDetectionThreshold,'r--')
        plot([0 length(data.vel)],[1 1]*data.saccadeVelocityTreshold,'r:')
        plot(data.sacon,data.vel(data.sacon),'bo','MarkerFaceColor','blue','MarkerSize',4)
        plot(data.sacoff,data.vel(data.sacoff),'ro','MarkerFaceColor','red','MarkerSize',4)
        qhighvelglissade = data.glissadetype==2;
        plot(data.glissadeoff(~qhighvelglissade),data.vel(data.glissadeoff(~qhighvelglissade)),'g*')
        plot(data.glissadeoff( qhighvelglissade),data.vel(data.glissadeoff( qhighvelglissade)),'c*')
        for p=1:length(data.sacoff)
            plot([0 ceil(ETparams.glissadeSearchWindowms/1000 * ETparams.samplingFreq)]+data.sacoff(p),[1 1]*data.localSaccadeVelocityTreshold(p),'r-');
        end
        linkaxes([ax gca]);
        % end temp

        % Implicitly detect fixations
        %-------------------------------------            
        detectFixations(i,j)

    end
end