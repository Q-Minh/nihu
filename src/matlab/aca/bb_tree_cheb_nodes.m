function xCheb = bb_tree_cheb_nodes(CTree, nExp, dim)
%BB_TREE_CHEB_NODES Chebyshev nodes of a cluster tree

nNode = nExp^dim;
nClusters = length(CTree);
xCheb = zeros(nClusters*nNode,dim);

X = repmat([-1;1], 1, dim);
x0 = chebroots(nExp, dim);

for c = 1 : nClusters
    idx = c*nNode + (1:nNode);
    xCheb(idx,:) = lintrans(x0, X, CTree(c).bb);
end

end
