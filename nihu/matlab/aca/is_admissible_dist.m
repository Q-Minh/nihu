function b = is_admissible_dist(C1, C2, eta)
bb1 = C1.bb;	% bounding box of cluster 1
bb2 = C2.bb;	% bounding box of cluster 2
d1 = norm(diff(bb1));	% diameter of cluster 1
d2 = norm(diff(bb2));	% diameter of cluster 2
% estimate minimal distance between two nodes
dist = max(norm(mean(bb1)-mean(bb2))-(d1+d2)/2, 0);
b = min(d1, d2) < eta * dist;
end % of function is_admissible
