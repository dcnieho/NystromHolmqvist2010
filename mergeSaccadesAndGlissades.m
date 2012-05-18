function [data] = mergeSaccadesAndGlissades(data, ETparams, extraCut)
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


% next, stretch up the part cut around saccade onsets and offsets
if nargin>2 && ~isempty(extraCut) && any(extraCut)
    data.saccade.on  = data.saccade.on  + ceil(extraCut(1)/1000 * ETparams.samplingFreq);
    data.saccade.off = data.saccade.off + ceil(extraCut(2)/1000 * ETparams.samplingFreq);
    
    % We could now have saccade starts before the end of the
    % previous saccade. prune/merge them here.
    data.saccade = mergeIntervals(data.saccade,[],0);
    
    % make sure first onset and last offset doesn't run out of the data
    if data.saccade.on(1) < 1
        data.saccade.on(1) = 1;
    end
    if data.saccade.off(end) > size(data.deg.vel,1)
        data.saccade.off(end) = size(data.deg.vel,1);
    end
end