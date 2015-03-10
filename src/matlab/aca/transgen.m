function trans = transgen(dim, n)
%TRANSGEN Generate unique translation vectors
%   TRANS = TRANSGEN(DIM, N) generates unique translation vectors in DIM
%   dimensions so that each coordinate can take n different values.s
%   The resultiung tarnslation vectors are symmetric to the origin.
%
% Example:
%    transgen(2, 3)
% results in the translation vectors
%     -1    -1
%      0    -1
%      1    -1
%     -1     0
%      0     0
%      1     0
%     -1     1
%      0     1
%      1     1
%
% See also: trans2idx

% Copyright (c) 2014 Peter Fiala

idx = (0 : n^dim-1)';
trans = zeros(length(idx),dim);
offset = (n-1)/2;
for d = 1 : dim
    trans(:,d) = mod(idx, n)-offset;
    idx = floor(idx / n);
end

end
