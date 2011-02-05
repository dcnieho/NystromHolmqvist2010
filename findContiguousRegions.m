function [on,off] = findContiguousRegions(in)
% finds all contiguous sections of true in a boolean array
% if a non-boolean is passed in, it is first converted to boolean, so
% anything nonzero will become true (excpet nan, which will error)
%
% Example:
% [on,off]=findContiguousRegions([1 0 0 0 1 1 1 0 1 1 0])
% on  =
%      1     5     9
% off =
%      1     7    10
%
% % you can then reconstruct the boolean vector with:
% idx = bounds2bool(on,off)
% idx =
%      1     0     0     0     1     1     1     0     1     1
%
% See also bounds2bool

if ~islogical(in)
    in = logical(in);
end


D   = diff([false in false]);
on  = find(D == 1);
off = find(D == -1) - 1;