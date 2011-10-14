function [mainSet,extraSet] = mergeIntervals(mainSet,extraSet,windowSamples)
% Merges intervals whose offset is closer than windowSamples to the next
% interval's onset (e.g. two saccades that are very close).
% If the second input is non-empty, it is checked if any onset in the extra
% set coincides with an offset in the main set. If so the corresponding
% offset of the extra set is used, instead of the offset of the main set,
% when checking the nearness of the next onset in the main set (e.g
% checking if the glissade offset is close to the onset of the next
% saccade).
% When windowSamples is set to 0, overlapping intervals are merged.

kk=1;
while kk < length(mainSet.on)      % NB: doesn't process last interval (useless anyway of course!)
    % walk through all saccades and see if followed shortly by another
    % saccade.
    
    % If there is a glissades for this saccade, check from the end of the
    % glissade, otherwise check from the end of the saccade
    if ~isempty(extraSet)
        qHaveExtraOffset = mainSet.off(kk)==extraSet.on;
    else
        qHaveExtraOffset = false;
    end
    if any(qHaveExtraOffset)
        assert(sum(qHaveExtraOffset)==1)  % any more would be ridiculous!
        thisoff = extraSet.off(qHaveExtraOffset);
    else
        thisoff = mainSet.off(kk);
    end
    
    % check if start of next saccade occurs within the window
    if mainSet.on(kk+1)-thisoff <= windowSamples
        % if yes, merge
        mainSet.off(kk) = mainSet.off(kk+1);
        
        % remove next saccade...
        mainSet         = replaceElementsInStruct(mainSet,kk+1,[]);
        
        % ... and glissade that is caught in between
        if ~isempty(extraSet)
            extraSet    = replaceElementsInStruct(extraSet,qHaveExtraOffset,[]);
        end
        
        continue;
    else
        % else, process next saccade
        kk = kk+1;
    end
end