function [S, b, nnz] = saistruct(nodes, d)
%SAISTRUCT  Determine sparsity structure of SAI preconditioner
%   [S, b, nnz] = saistruct(nodes, d) determines the block diagonal
%   sparsity structure of a preconditioner matrix. The input nodes are
%   clustered using a cluster diameter d, and full submatrices are
%   associated to each cluster
% Parameters:
%   nodes : Nx3 xyz coordinates of nodes
%   d     : cluster diameter
%   S     : Bx? Sparsity structure containing B blocks. The matrix contains
%           the node indices of a cluster in each row.
%   b     : Bx1 block sizes
%   nnz   : number of nonzero entries (dot(b,b))
%
% See also: saileft, sairight

% Peter Fiala
% 2009

%%
relnodes = floor((nodes - repmat(min(nodes,[],1), size(nodes,1), 1)) / d);
[relnodes, m, n] = unique(relnodes, 'rows');
S = sort(sparse(n, 1:length(n), 1:length(n)), 2);
S = fliplr(full(S(:,any(S,1))));
b = sum(S ~= 0, 2);
nnz = dot(b, b);
end