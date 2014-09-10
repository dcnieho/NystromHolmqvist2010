function [data] = mergeSaccadesAndGlissades(data, ETparams)
% NB: You would need to rerun processSaccadesAndGlissades.m after this, if
% you are interested in its output, as the saccades flags have just changed
%
% Also, this function does not remove the glissade markers from the data
% struct, you have to do that yourself if you no longer need them:
% if isfield(data,'glissade')
%     data = rmfield(data,'glissade');    % remove glissades, we fused them anyway
% end

% get saccade offset markers
sacoff  = data.saccade.off;

if isempty(sacoff)
    % nothing to do
    return
end

if isfield(data,'glissade')
    for p = 1:length(data.glissade.off)
        % for each glissade, find corresponding saccade
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
