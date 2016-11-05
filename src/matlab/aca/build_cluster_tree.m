function Ctree = build_cluster_tree(x, Cleaf, divtype)
%BUILD_CLUSTER_TREE build cluster tree from a set of coordinates
%   T = BUILD_CLUSTER_TREE(x) builds a cluster tree containing a
%   hierarchical structure of clusters. The first cluster encapsulates
%   each point, and the leaf clusters contain a limited number of points.
%
%   T = BUILD_CLUSTER_TREE(X, C) creates a cluster tree where the leafs
%   contain max C nodes. The default for C is 10.
%
%   T = BUILD_CLUSTER_TREE(X, C, DIVTYPE) creates a cluster tree where the leafs
%   contain max C nodes. The default for C is 10.
%
% See also: build_block_tree
%
% Example:
%   m = create_sphere_boundary(1, 10);
%   C = build_cluster_tree(m.Nodes(:,2:4));
%
% Copyright (C) 2016 Peter Fiala

if nargin < 3
    divtype = 'svd';
end

if nargin < 2
    Cleaf = 10;
end

switch divtype
    case 'svd'
        divide_fun = @(C,x)get_svd_child_clusters(C, x);
    case 'oc'
        divide_fun = @(C,x)get_oc_child_clusters(C, x);
    otherwise
        error('build_cluster_tree:invalid_argument', ...
            'Argument ''divtype'' must be ''svd'' or ''oc''');
end

bb = [min(x,[],1); max(x,[],1)];
if strcmp(divtype, 'oc')
    c = mean(bb,1);
    D = max(diff(bb,1)) * (1+1e-2);
    bb = bsxfun(@plus, c, [-D; D]/2 * ones(1,size(x,2)));
end

Ctree(1).ind = 1 : size(x,1);
Ctree(1).level = 0;
Ctree(1).children = [];
Ctree(1).child_idx = [];
Ctree(1).bb = bb;

Capacity = 100;
Ctree(Capacity).level = 0;
i = 1;
S = 1;
while i <= S
    if length(Ctree(i).ind) > Cleaf
        [Children, Idx] = divide_fun(Ctree(i), x);
        idx = S + (1 : length(Children));
        Ctree(i).children = idx;
        Ctree(i).child_idx = Idx;
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

for i = 1 : length(Ctree)
    for j = 1 : length(Ctree(i).children)
        Ctree(Ctree(i).children(j)).father = i;
        Ctree(Ctree(i).children(j)).father_idx = Ctree(i).child_idx(j);
    end
end

end

function [CC, Idx] = get_svd_child_clusters(C, x)

x0 = x(C.ind,:);
[~,~,V] = svd(x0, 0);

s = std(x0 * V);
[~, i] = max(s);
dir = V(:,i);

[~, i] = sort(x0 * dir);
nhalf = round(length(i)/2);

CC(1).ind = C.ind(i(1:nhalf));
CC(2).ind = C.ind(i(nhalf+1:end));

for k = 1 : 2
    CC(k).level = C.level+1;
    CC(k).children = [];
    CC(k).bb = [min(x(CC(k).ind,:)); max(x(CC(k).ind,:))];
end

end

function [CC, Idx] = get_oc_child_clusters(C, x)
dim = size(x,2);
c0 = mean(C.bb);
D = diff(C.bb,1);
D = D(1);
Idx = (0 : 2^dim-1);
d = ((dec2bin(Idx)-'0')-.5) * D/2;
n = size(d,1);
c = bsxfun(@plus, c0, d);
CC = struct(...
    'ind', mat2cell(repmat(C.ind, n, 1), ones(n,1), length(C.ind)), ...
    'level', num2cell(repmat(C.level+1, n, 1)), ...
    'children', mat2cell(zeros(n,0), ones(n,1), 0), ...
    'child_idx', mat2cell(zeros(n,0), ones(n,1), 0));

include = true(n,1);
for i = 1 : length(CC)
    CC(i).bb = bsxfun(@plus, c(i,:), [-D/4; D/4] * ones(1,dim));
    CC(i).ind = CC(i).ind(all(...
        bsxfun(@gt, x(CC(i).ind,:), CC(i).bb(1,:)) &...
        bsxfun(@le, x(CC(i).ind,:), CC(i).bb(2,:)),...
        2));
    include(i) = ~isempty(CC(i).ind);
end
CC = CC(include);
Idx = Idx(include);
end
