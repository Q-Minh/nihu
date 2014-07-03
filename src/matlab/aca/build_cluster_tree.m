function Ctree = build_cluster_tree(x, Cleaf)
%BUILD_CLUSTER_TREE build cluster tree for ACA algorithm
%   T = BUILD_CLUSTER_TREE(x) builds a cluster tree containing a
%   hierarchical structure of clusters. The first cluster encapsulates
%   each point, and the leaf clusters contain one single point only.
%
%   T = BUILD_CLUSTER_TREE(X, C) creates a cluster tree where the leafs
%   contain max C nodes. The default for C is 10.
%
% See also: build_block_tree
%
% Example:
%   m = create_sphere_boundary(1, 10);
%   C = build_cluster_tree(m.Nodes(:,2:4));
%
% Copyright (C) 2014 Peter Fiala

if nargin == 1
    Cleaf = 10;
end

Ctree(1).ind = 1 : size(x,1);
Ctree(1).level = 0;
Ctree(1).children = [];

Capacity = 100;
Ctree(Capacity).level = 0;
i = 1;
S = 1;
while i <= length(Ctree)
    if length(Ctree(i).ind) <= Cleaf
    else
        Children = get_child_clusters(Ctree(i), x);
        idx = S + (1 : length(Children));
        Ctree(i).children = idx;
        if (max(idx) > Capacity)
            Capacity = 2*Capacity;
            Ctree(Capacity).level = 0;
        end
        Ctree(idx) = Children;
        S = S + length(Children);
    end
    i = i + 1;
end
Ctree = Ctree(1:S);

end

function CC = get_child_clusters(C, x)

x0 = x(C.ind,:);
[~,~,V] = svd(x0, 0);

s = std(x0 * V);
[~, i] = max(s);
dir = V(:,i);

[~, i] = sort(x0 * dir);
nhalf = round(length(i)/2);

CC(1).ind = C.ind(i(1:nhalf));
CC(1).level = C.level+1;
CC(1).children = [];

CC(2).ind = C.ind(i(nhalf+1:end));
CC(2).level = C.level+1;
CC(2).children = [];

end
