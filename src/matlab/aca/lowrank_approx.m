function [U, V] = lowrank_approx(M, siz, eps, R)
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

rows = siz(1);
cols = siz(2);

if nargin < 4
    R = min([rows cols]);
end

U = zeros(rows,R);
V = zeros(cols,R);
S2 = 0;

indices = [];

i = 1;
r = 0;
while r < R
    row = M(i,1:cols) - U(i,1:r) * V(:,1:r)';
    [gamma, j] = max(abs(row));
    if gamma < 1e-10
        i = setdiff(1:rows, indices);
        if (isempty(i))
            U = U(:, 1:r);
            V = V(:, 1:r);
            return;
        end
        i = i(1);
        continue;
    end
    
    indices(end+1) = i;
    
    col = M(1:rows,j) - U(:,1:r) * V(j,1:r)';
    
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
    
    [~,i] = max(abs(col(setdiff(1:rows,indices))));
end

end
