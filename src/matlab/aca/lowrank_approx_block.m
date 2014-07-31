function [U, V] = lowrank_approx_block(M, I, J, eps)
%LOWRANK_APPROX return low rank approximation of a matrix block
%   [U,V] = LOWRANK_APPROX(M, I, J, eps) returns the low-rank approximation
%   of the matrix block M(I,J). The low rank approximation is of the form
%   M(I:J) = U * V
%   where U is an m x R matrix and V is an R x n matrix.
%   Argument R denotes the maximal rank of the approximation.
%
% See also: BUILD_CLUSTER_TREE, BUILD_BLOCK_TREE
%
% Copyright (C) 2014 Peter Fiala

% block
B = @(i,j)M(I(i), J(j));

[U, V] = lowrank_approx(B, [length(I), length(J)], eps);

end
