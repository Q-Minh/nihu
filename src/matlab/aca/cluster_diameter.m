function D = cluster_diameter(Ctree, x)
nTree = length(Ctree);
D = zeros(nTree,1);

for c = 1 : nTree
    x0 = x(Ctree(c).ind,:);
    [~, ~, V] = svd(x0, 0);
    x0aligned = x0 * V;
    D(c) = norm(max(x0aligned)-min(x0aligned));
end
end
