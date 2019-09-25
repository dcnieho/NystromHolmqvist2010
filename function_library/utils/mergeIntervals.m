function [mainSet,extraSet] = mergeIntervals(mainSet,extraSet,windowSamples,qNaN)
% Merges intervals whose offset is closer than windowSamples to the next
% interval's onset (e.g. two saccades that are very close).
% If the second input is non-empty, it is checked if any onset in the extra
% set coincides with an offset in the main set. If so the corresponding
% offset of the extra set is used, instead of the offset of the main set,
% when checking the nearness of the next onset in the main set (e.g
% checking if the glissade offset is close to the onset of the next
% saccade).
% When windowSamples is set to 0, only overlapping intervals are merged.
% the code below is not commented in general words, but for the case of
% merging saccades, which might be followed by glissades as that's easier
% to follow.

if nargin<4
    qDontCountNan = false;
else
    qDontCountNan = true;
end

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
    if (  qDontCountNan && sum(~qNaN(thisoff:mainSet.on(kk+1))) <= windowSamples) || ...    % check less than windowSamples non-nan intervening samples
        (~qDontCountNan && mainSet.on(kk+1)-thisoff <= windowSamples)                       % or just check number of  intervening samples
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