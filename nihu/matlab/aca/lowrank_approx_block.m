function [U, V] = lowrank_approx_block(M, I, J, varargin)
%LOWRANK_APPROX_BLOCK low rank approximation of a matrix block
%   [U,V] = LOWRANK_APPROX(M, I, J, eps, R) returns the low-rank
%   approximation of the matrix block M(I,J).
%   The low rank approximation is of the form M(I,J) = U * V'
%   Argument R denotes the maximal rank of the approximation.
%
% See also: LOWRANK_APPROX
%
% Copyright (C) 2014 Peter Fiala

[U, V] = lowrank_approx(@(i,j)M(I(i), J(j)), [length(I), length(J)],...
    varargin{:});
end
