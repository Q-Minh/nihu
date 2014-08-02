function [clusters, perm, order] = sort_cluster_tree(CTree)
%SORT_CLUSTER_TREE cluster-contiguous node numbering
%   [C, P, Q] = sort_cluster_tree(Tree)
%   generates cluster-contiguous node numbering for a cluster tree.
%   If the tree was generated from the nodal locations X, then
%   X(P,:) is a cluster-contiguous permutation of the locations.
%
%   Each row of matrix C contains the offset (starting with 0) and size
%   of a contiguous cluster.
%
%   Q is the inverse permutation of P.
%
% see also: build_cluster_tree, build_block_tree
%
% Example:
%  xs = rand(1000,2);
%  T = build_cluster_tree(xs,100);
%  [C, P] = sort_cluster_tree(T);
%  range = C(end,1) + (1:C(end,2));     % last cluster
%  plot(xs(:,1), xs(:,2), 'k.', xs(P(range),1), xs(P(range),2), 'r*');

% Copyright (C) 2014 Peter Fiala

clusters = zeros(length(CTree),2);

% get total number of indices
isleaf = false(length(CTree),1);
for t = 1 : length(CTree)
    if isempty(CTree(t).children)
        isleaf(t) = true;
    end
    clusters(t,2) = length(CTree(t).ind);
end

leafs = CTree(isleaf);
perm = [];
for t = 1 : length(leafs)
    perm = [perm leafs(t).ind];
end
order(perm) = 1 : length(perm);

for t = 1 : length(CTree)
    clusters(t,1) = min(order(CTree(t).ind))-1;
end

end
