function [U, V] = lowrank_approx(M, eps, R)
%LOWRANK_APPROX return low rank approximation of a matrix
%   [U,V] = LOWRANK_APPROX(M, eps) returns the low-rank approximation
%   of the matrix block M. The low rank approximation is of the form
%   M = U * V
%   where U is an m x R matrix and V is an R x n matrix.
%   Argument R denotes the maximal rank of the approximation.
%
% See also: BUILD_CLUSTER_TREE, BUILD_BLOCK_TREE
%
% Copyright (C) 2014 Peter Fiala

[m, n] = size(M);

if nargin < 3
    R = min(m,n);
end

U = zeros(m,R);
V = zeros(R,n);

i = 0;
r = 0;
while true
    i = i + 1;
    if i > m
        break;
    end
    
    row = M(i,:) - U(i,1:r) * V(1:r,:);
    [g, j] = max(abs(row));
    if g < 1e-8
        continue;
    end
    col = M(:,j) - U(:,1:r) * V(1:r,j);
    
    r = r + 1;
    
    U(:, r) = col/col(i);
    V(r, :) = row;
    
    if norm(U(:,r), 'fro') * norm(V(r,:), 'fro') < eps * norm(U*V, 'fro')
        break;
    end
end

U = U(:, 1:r);
V = V(1:r, :);

end
