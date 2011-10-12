function [data] = mergeSaccadesAndGlissades(data)

% preprocess on and offsets:
% - fuse glissade after saccade
% - the two above also incorporates: fuse saccade that starts at glissade offset

% get saccade onset and offset markers
sacon   = data.saccade.on;
sacoff  = data.saccade.off;

if isfield(data,'glissade')
    for p = 1:length(data.glissade.off)
        % merge all glissades first before merging saccades, otherwise its
        % difficult to handle complicated cases where two saccades, both
        % followed with glissades, need to be fused (or even nastier).
        % Splitting up is easier and robust.
        
        % for eah glissade, find corresponding saccade
        % easy, as glissade onset is equal to fixation offset per definition
        qsac = data.glissade.on(p) == sacoff;
        assert(sum(qsac)==1)
        
        % change saccade offset to glissade offset
        sacoff(qsac) = data.glissade.off(p);
    end
end

data.saccade.on     = sacon;
data.saccade.off    = sacoff;


% kk=1;
% while kk < length(sacon)    % NB: doesn't process last saccade (useless of course!)
%     % walk through all saccades and see if followed shortly by another
%     % saccade
%     
%     if sacon(kk+1)-sacoff(kk) <= SacMergeWindowSamp
%         % if yes, merge, continue
%         sacoff(kk) = sacoff(kk+1);
%         
%         sacon(kk+1)  = [];
%         sacoff(kk+1) = [];
%         
%         continue;
%     else
%         % else, process next saccade
%         kk = kk+1;
%     end
% end