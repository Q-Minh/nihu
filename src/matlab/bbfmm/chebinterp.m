function S = chebinterp(n, y, Y)
%CHEBINTERP Chebyshev anterpolation matrix
%    S = CHEBINTERP(n, y, Y) generates the anterpolation matrix that
%    anterpolates from points y to n Chebyshev nodes in interval Y.
%
% Example:
%    S = chebinterp(5, rand(100,2), [0 0; 2 2]);
% anterpolates from 100 random 2D locations to 5x5 Chebyshev nodes in the
% [0;2] x [0;2] interval

narginchk(2,3);

d = size(y,2);  % dimensionality

if nargin == 2
    Y = repmat([-1 1]', 1, d);
end

if size(Y,2) ~= d
    error('bbFMM:arg', 'bounding box and input dimensions must match');
end

x = chebroots(n, d);

y = lintrans(y, Y, repmat([-1 1]', 1, d));
S = ones(size(x,1), size(y,1));
for j = 1 : d
    S = S .* chebinterp0(n, x(:,j), y(:,j));
end

end % of function CHEBINTERP

function S = chebinterp0(N, x, y)
n = length(x);
ny = length(y);
S = ones(n,ny);
for l = 1 : N-1
    S = S + 2*chebpoly(l, x) * chebpoly(l, y');
end
S = S / N;
end
