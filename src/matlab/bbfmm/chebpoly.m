function T = chebpoly(n, x)
%CHEBPOLY values of Chebyshev polynomials
%  T = CHEBPOLY(N, X) computes the values of the Chebyshev polynomial Tn(x)
%  of order N at the locations X
%  N and X should be of the same size
% Example:
%   n = (1:10)';
%   x = -1 : 1e-2 : 1;
%   plot(x, chebpoly(n, x));
T = cos(bsxfun(@times, n, acos(x)));
end
