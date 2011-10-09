function [output] = rovingBin(input,bin_width,pre_post)

% function [output] = rovingBin(input,bin_width,pre_post);
%
% this function outputs an array containing length(input)
% rows and bin_width columns
%
% pre_post specifies how many columns there are before and after
% each entry.

% roving_bin(1:10,10,[9 0])
%
% ans =
%
%    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1
%    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     2
%    NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     2     3
%    NaN   NaN   NaN   NaN   NaN   NaN     1     2     3     4
%    NaN   NaN   NaN   NaN   NaN     1     2     3     4     5
%    NaN   NaN   NaN   NaN     1     2     3     4     5     6
%    NaN   NaN   NaN     1     2     3     4     5     6     7
%    NaN   NaN     1     2     3     4     5     6     7     8
%    NaN     1     2     3     4     5     6     7     8     9
%      1     2     3     4     5     6     7     8     9    10


if size(input,1)~=1
    input = input';
end

if nargin<3
    
    if mod(bin_width,2) ==1 % odd number means each bin is centered
        
        pre_bin_width = floor(bin_width/2);
        post_bin_width = floor(bin_width/2);
        
        %pre_bin_width  = bin_width -1;
        %post_bin_width = 0;
        
    elseif mod(bin_width,2) == 0 % even bin, 'off-center' the value behind
        
        pre_bin_width  = bin_width/2;
        post_bin_width = bin_width/2 - 1;
        
    end
    
else
    assert(sum(pre_post)+1==bin_width,'pre_post length not compatible with bin width')
    pre_bin_width  = pre_post(1);
    post_bin_width = pre_post(2);
end

if 0
    output = zeros(length(input),bin_width);
    
    % center the observations within the bin
    
    pre_bin = zeros(1,pre_bin_width);
    pre_bin(pre_bin==0) = nan;
    
    post_bin = zeros(1,post_bin_width);
    post_bin(post_bin==0) = nan;
    
    padded_input = [pre_bin input post_bin];
    
    % if this could be done without a loop,that would be better
    
    for i=pre_bin_width+1:pre_bin_width+length(input)
        % start at the first data entry
        % through the last point of data
        
        output(i-pre_bin_width,:) = padded_input(i-pre_bin_width:i+post_bin_width);
        
    end
else
    % prepare output
    output = nan(length(input),bin_width);
    
    % make index vector
    x = [-pre_bin_width : post_bin_width];
    ind = repmat([1:length(input)].',1,bin_width) + repmat(x,length(input),1);
    
    % ignore indices falling outside of data
    q = (ind>0 & ind<=length(input));
    
    % make data matrix
    output(q) = input(ind(q));
end