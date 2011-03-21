function [sacon,sacoff] = MergeSaccadesAndGlissades(data,ETparams,SacMergeWindow)

% preprocess on and offsets:
% - fuse glissade after saccade
% - fuse two saccades with little time between them
% - the two above also incorporates: fuse saccade that starts at glissade offset

% if saccade is followed by another saccade by less than SacMergeWindow (in
% ms), they'll be merged
SacMergeWindowSamp  = ceil(SacMergeWindow./1000 * ETparams.samplingFreq);

% get saccade onset and offset markers
sacon = data.saccade.on;
sacoff = data.saccade.off;
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


kk=1;
while kk < length(sacon)    % NB: don't process last saccade (useless of course!)
    % walk through all saccades and see if followed shortly by another
    % saccade
    
    if kk==25
        3;
    end
    
    if sacon(kk+1)-sacoff(kk) <= SacMergeWindowSamp
        % if yes, merge, continue
        sacoff(kk) = sacoff(kk+1);
        
        sacon(kk+1)  = [];
        sacoff(kk+1) = [];
        
        continue;
    else
        % else, process next saccade
        kk = kk+1;
    end
end