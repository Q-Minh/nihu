function Ctree = build_oc_cluster_tree(x, Cleaf)
%BUILD_OC_CLUSTER_TREE build cluster octree for ACA algorithm
%   T = BUILD_CLUSTER_TREE(x) builds a cluster tree containing a
%   hierarchical structure of clusters. The first cluster encapsulates
%   each point, and the leaf clusters contain one single point only.
%
% See also: build_block_tree
%
% Example:
%   m = create_sphere_boundary(1, 10);
%   C = build_oc_cluster_tree(m.Nodes(:,2:4));
%
% Copyright (C) 2014 Peter Fiala

if nargin == 1
    Cleaf = 10;
end

bb = [min(x,[],1); max(x,[],1)];

Ctree(1).ind = 1 : size(x,1);
Ctree(1).level = 0;
Ctree(1).children = [];
Ctree(1).c = mean(bb, 1);
Ctree(1).D = max(diff(bb, 1)) + 1e-3;

Capacity = 100;
Ctree(Capacity).level = 0;
i = 1;
S = 1;
while i <= S
    if length(Ctree(i).ind) > Cleaf
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
d = [
    -1 -1 -1
    1 -1 -1
    1  1 -1
    -1  1 -1
    -1 -1  1
    1 -1  1
    1  1  1
    -1  1  1
    ] * C.D/4;
c = bsxfun(@plus, C.c, d);
CC = struct(...
    'ind', mat2cell(repmat(C.ind, 8, 1), ones(8,1), length(C.ind)), ...
    'level', num2cell(repmat(C.level+1, 8, 1)), ...
    'children', mat2cell(zeros(8,0), ones(8,1), 0),...
    'c', mat2cell(c, ones(8,1), 3),...
    'D', C.D/2);

include = true(8,1);
for i = 1 : length(CC)
    CC(i).ind = CC(i).ind(all(bsxfun(@gt, x(CC(i).ind,:), CC(i).c-CC(i).D/2) &...
        bsxfun(@le, x(CC(i).ind,:), CC(i).c+CC(i).D/2), 2));
    include(i) = ~isempty(CC(i).ind);
end
CC = CC(include);
end
