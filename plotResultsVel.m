function plotResultsVel(ETparams,i_idx,j_idx)
%--------------------------------------------------------------------------
% Process eye-movement data for file and trial separately, i - files, j -
% trials
%--------------------------------------------------------------------------
if nargin < 2
    i_idx = 1:size(ETparams.data,1);
    j_idx = 1:size(ETparams.data,2);
end

h1 = figure;
for i = i_idx
    for j = j_idx
        
        clf
        hold on
        
        % h1 = plotTrialDetectionResults(ETparams,i,j,h1,'r',-5);
        h1 = plotTrialVelDetection(ETparams,i,j,h1,'r',-2);
        hold off
%         pause
        
    end
end

