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
V = zeros(n,R);
S2 = 0;

i = 1;
r = 0;
while r < R
    row = M(i,:) - U(i,1:r) * V(:,1:r)';
    [gamma, j] = max(abs(row));
    if gamma < 1e-8
        i = i + 1;
        continue;
    end
    col = M(:,j) - U(:,1:r) * V(j,1:r)';
    
    r = r + 1;
    
    U(:, r) = col/col(i);
    V(:, r) = row';
    
    S2 = S2 + norm(U(:,r), 'fro')^2 * norm(V(:,r), 'fro')^2;
    for j = 1 : r-1
        S2 = S2 + 2 * (U(:,r)'*U(:,j)) * (V(:,j)'*V(:,r));
    end
    
    if norm(U(:,r), 'fro') * norm(V(:,r), 'fro') < eps * sqrt(S2)
        U = U(:, 1:r);
        V = V(:, 1:r);
        return;
    end
    
    [~,i] = max(abs(col(setdiff(1:m,i))));
end

end
