function [U, V] = lowrank_approx(M, I, J, eps)
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

m = length(I);
n = length(J);

% if m == 1
%     U = 1;
%     V = M(I,J);
%     return;
% end
% 
% if n == 1
%     U = M(I,J);
%     V = 1;
%     return;
% end
% 
R = 20;

U = zeros(m,R);
V = zeros(R,n);

i = 0;
r = 0;
while true
    i = i + 1;
    if i > m
        break;
    end
    
    row = M(I(i),J);
    row = row - U(i,1:r) * V(1:r,:);
    [g, j] = max(abs(row));
    if g < 1e-8
        continue;
    end
    col = M(I,J(j));
    col = col - U(:,1:r) * V(1:r,j);
    
    r = r + 1;
    
    U(:, r) = col/col(i);
    V(r, :) = row;
    
    if norm(U(:,r)*V(r,:), 'fro') / norm(U*V, 'fro') < eps
        break;
    end
end

U = U(:, 1:r);
V = V(1:r, :);

end
