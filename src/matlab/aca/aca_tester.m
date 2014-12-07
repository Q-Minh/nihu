function [times, res] = aca_tester(M, xr, xs, exc, nLeaf, Admit, eps)
%ACA_TESTER general tester of the ACA algorithm

Ns = size(xs,1);
Nr = size(xr,1);

t0 = tic;
SourceTree = build_cluster_tree(xs, nLeaf);
ReceiverTree = build_cluster_tree(xr, nLeaf);
times.tTree = toc(t0);

t0 = tic;
[B_near, B_far] = build_dual_block_tree(ReceiverTree, SourceTree, Admit);
times.tBtree = toc(t0);

res.resp = zeros(Nr, 1);

t0 = tic;
for b = 1 : size(B_near,1)
    I = ReceiverTree(B_near(b,1)).ind;
    J = SourceTree(B_near(b,2)).ind;
    res.resp(I) = res.resp(I) + M(I,J) * exc(J);
end
times.tNear = toc(t0);

t0 = tic;
for b = 1 : size(B_far,1)
    progbar(1, size(B_far,1), b);
    I = ReceiverTree(B_far(b,1)).ind;
    J = SourceTree(B_far(b,2)).ind;
    [U, V] = lowrank_approx_block(M, I, J, eps);
    res.resp(I) = res.resp(I) + U * (V' * exc(J));
end
times.tAca = toc(t0);

t0 = tic;
res.resp0 = M(1:Nr,1:Ns) * exc;
times.tFull = toc(t0);

end
