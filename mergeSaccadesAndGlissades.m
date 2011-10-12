function [data] = mergeSaccadesAndGlissades(data)

% get saccade offset markers
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

% store merged offset markers. Don't remove glissade markers, user can do
% that themselves if they don't want them.
data.saccade.off    = sacoff;