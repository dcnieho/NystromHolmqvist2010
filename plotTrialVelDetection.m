function h1 = plotTrialVelDetection(ET,i,j,h1,c1,offset)
% A velocity plot with detection results indicated


convf = 1/(ET.samplingFreq/1000); % Conversion factor between samples and ms

%--------------------------------
% Parameters 
%--------------------------------
X = ET.data(i,j).X;
Y = ET.data(i,j).Y;
V = ET.data(i,j).vel;
A = ET.data(i,j).acc;

nSamples = length(V);

%--------------------------------
% Plot XY-velocities
%--------------------------------   
if strcmp(c1,'r') == 1
    plot((1:length(V))*convf,V,c1,'LineWidth',2);
else
    plot((1:length(V))*convf,V,c1,'LineWidth',2);
end

%--------------------------------
% Indicate periods of fixation and saccade
%--------------------------------   

% Plot fixations
idx = find(ET.fixationIdx(i,j).Idx);
plot(idx*convf,offset*ones(1,length(idx)),[c1,'.'],'LineWidth',5);

% Plot detection thresholds
plot((1:nSamples)*convf,ones(1,nSamples)*ET.data(i,j).saccadeVelocityTreshold,'-'...
    ,(1:nSamples)*convf,ones(1,nSamples)*ET.data(i,j).peakDetectionThreshold,':');

Vmax = 500;
% Plot saccades
for kk = 1:size(ET.saccadeInfo,3)
    idxSaccadeEnd = ET.saccadeInfo(i,j,kk).end*ET.samplingFreq;
    idxSaccadeStart = ET.saccadeInfo(i,j,kk).start*ET.samplingFreq;
    
    if ~isempty(idxSaccadeEnd )
        plot(ones(1,length(0:Vmax))*idxSaccadeStart*convf,0:Vmax,'-','LineWidth',1);
        plot(ones(1,length(0:Vmax))*idxSaccadeEnd*convf,0:Vmax,'--','LineWidth',1);
        % Plot local threshold
        plot((idxSaccadeStart:idxSaccadeEnd)*convf,ones(1,length(idxSaccadeStart:idxSaccadeEnd))*ET.data(i,j,kk).localSaccadeVelocityTreshold,'-k','LineWidth',1);        
    end   
end

% Plot glissades
for kk = 1:size(ET.glissadeInfo,3)
    idxGlissadeEnd = ET.saccadeInfo(i,j,kk).end*ET.samplingFreq +  ET.glissadeInfo(i,j,kk).duration*ET.samplingFreq;
    if ~isempty(idxGlissadeEnd )
        plot(ones(1,length(0:Vmax))*idxGlissadeEnd*convf,0:Vmax,'-.','LineWidth',1);
    end
end
title(['Subject: ',num2str(i),' Trial: ',num2str(j)])


axis([0 nSamples -20 max(V)]), axis xy
xlabel('Time (ms)')
ylabel('Velocity (${}^{\circ/s}$)','Interpreter','latex')



