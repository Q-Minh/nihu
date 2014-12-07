function p = display_block_structure(RowClusters, ColClusters, BlockTree)
%DISPLAY_BLOCK_STRUCTURE visualise the block tree
%  DISPLAY_BLOCK_STRUCTURE(RC, CC, B) plots the block tree structure
%  defined by the row clusters RC, column clusters CC and block tree B.
%  P = DISPLAY_BLOCK_STRUCTURE(...) returns the handle to the plotted patch
%  object.
%
% Example:
%  source = create_line(1, 5000);
%  receiver = translate_mesh(create_line(1, 5000), [.5 0 0]);
%  
%  ColTree = build_cluster_tree(centnorm(source), 20);
%  CC = sort_cluster_tree(ColTree);
%  RowTree = build_cluster_tree(centnorm(receiver), 20);
%  RC = sort_cluster_tree(RowTree);
%  [~, B] = build_block_tree_2(RowTree, ColTree, .8);
%  
%  display_block_structure(RC, CC, B);
%
% see also: build_cluster_tree build_block_tree_2 sort_cluster_tree

% Copyright (c) Peter Fiala 2014

p = zeros(size(BlockTree,1),1);
for b = 1 : size(BlockTree,1)
    i = RowClusters(BlockTree(b,1),:);
    j = ColClusters(BlockTree(b,2),:);
	X = j(1) + j(2)*[0 1 1 0];
	Y = i(1) + i(2)*[0 0 1 1];
    p(b) = patch(X,Y,b*ones(size(X)));
end

set(gca, 'yDir', 'reverse');

end % of function display_block_structure
