function Ctree = build_cluster_tree(x, Cleaf)
%BUILD_CLUSTER_TREE build cluster tree for ACA algorithm
%   T = BUILD_CLUSTER_TREE(x) builds a cluster tree containing a
%   hierarchical structure of clusters. The first cluster encapsulates
%   each point, and the leaf clusters contain one single point only.
%
%   T = BUILD_CLUSTER_TREE(X, C) creates a cluster tree where the leafs
%   contain max C nodes. The default for C is 20.
%
% See also: build_block_tree
%
% Example:
%   m = create_sphere_boundary(1, 10);
%   C = build_cluster_tree(m.Nodes(:,2:4));
%
% Copyright (C) 2014 Peter Fiala

if nargin == 1
    Cleaf = 20;
end

Ctree(1).ind = 1 : size(x,1);
Ctree(1).bb = [
    min(x(Ctree(1).ind,:),[],1)
    max(x(Ctree(1).ind,:),[],1)
    ];
Ctree(1).D = sqrt(sum(diff(Ctree(1).bb, 1).^2));
Ctree(1).children = [];

i = 1;
while i <= length(Ctree)
    if length(Ctree(i).ind) <= Cleaf
    else
        Children = get_child_clusters(Ctree(i), x);
        idx = length(Ctree) + (1 : length(Children));
        Ctree(i).children = idx;
        Ctree(idx) = Children;
    end
    i = i + 1;
end

end

function CC = get_child_clusters(C, x)
[Dmax, maxdir] = max(diff(C.bb, 1));
sep = C.bb(1,maxdir)+Dmax/2;

CC(1) = C;
CC(1).ind = C.ind(x(C.ind,maxdir) < sep);
CC(2) = C;
CC(2).ind = setdiff(C.ind, CC(1).ind);

for i = 1 : 2
    CC(i).bb = [
        min(x(CC(i).ind,:),[],1)
        max(x(CC(i).ind,:),[],1)
        ];
    CC(i).D = sqrt(sum(diff(CC(i).bb, 1).^2));
    CC(i).children = [];
end
end
