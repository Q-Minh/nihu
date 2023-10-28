function [times, res] = aca_tester_mex(k, xr, xs, exc, nLeaf, Admit, eps)
%ACA_TESTER general tester of the ACA algorithm

M = @(i,j)helmholtz_kernel(xr(i,:),xs(j,:),k);

Ns = size(xs,1);    % number of source points
Nr = size(xr,1);    % number of receiver points

% build cluster trees
t0 = tic;
SourceTree = build_cluster_tree(xs, nLeaf);
[CS, PS] = sort_cluster_tree(SourceTree);
ReceiverTree = build_cluster_tree(xr, nLeaf);
[CR, PR] = sort_cluster_tree(ReceiverTree);
times.tTree = toc(t0);

% build blocks
t0 = tic;
[B_near, B_far] = build_dual_block_tree(ReceiverTree, SourceTree, Admit);
times.tBtree = toc(t0);

res.resp_near = zeros(Nr, 1);
res.resp_far = zeros(Nr, 1);
res.resp_far_mex = zeros(Nr, 1);

% direct evaluation of near contributions
t0 = tic;
for b = 1 : size(B_near,1)
    I = ReceiverTree(B_near(b,1)).ind;
    J = SourceTree(B_near(b,2)).ind;
    res.resp_near(I) = res.resp_near(I) + M(I,J) * exc(J);
end
times.tNear = toc(t0);

% aca approximation of far blocks
t0 = tic;
for b = 1 : size(B_far,1)
    progbar(1, size(B_far,1), b);
    I = ReceiverTree(B_far(b,1)).ind;
    J = SourceTree(B_far(b,2)).ind;
    [U, V] = lowrank_approx_block(M, I, J, eps);
    if b == size(B_far,1)
        disp(U);
        disp(V);
    end
    res.resp_far(I) = res.resp_far(I) + U * (V' * exc(J));
end
times.tAca = toc(t0);

% aca approximation of far blocks with mex
t0 = tic();
[res.resp_far_mex(PR), outranks] = aca_test(ensure_complex(k),...
    xs(PS,:), xr(PR,:),...
    int32(CS), int32(CR),...
    int32(B_far-1),...
    eps, int32(1000),...
    ensure_complex(exc(PS)));
times.tAcaMex = toc(t0);

% direct evaluation for comparison
t0 = tic;
res.resp0 = M(1:Nr,1:Ns) * exc;
times.tFull = toc(t0);

end
