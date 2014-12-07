function comp = max_aca_compression(rc, cc, bfar)
%maximal compression of the aca algorithm
%   COMP = MAX_ACA_COMPRESSION(RC, CC, BFAR) returns the maximal data
%   compression of the ACA algorithm based on the block structure BFAR, and
%   the row and column cluster trees RC and CC. The maximal data
%   compression is estimated by assuming that each block can be represented
%   by a diad.
%
% Example:
%   Ctree = build_cluster_tree(x, nLeaf);
%   [~, B_far] = build_block_tree(Ctree, Admit);
%   C = sort_cluster_tree(Ctree);
%   comp = max_aca_compression(C, C, B_far);

nTotal = dot(rc(bfar(:,1),2), cc(bfar(:,2),2), 1);
nMin = sum(2 * (rc(bfar(:,1),2) + rc(bfar(:,2),2)));
comp = nTotal/nMin;

end
