clear;

%% parameters
L = 1;
N = 5000;
nLeaf = 25;

%%
source = create_line(L, N);
receiver = translate_mesh(create_line(L, N), [L/2 0 0]);
xs = centnorm(source);
xr = centnorm(receiver);

%% build cluster trees
SourceTree = build_cluster_tree(xs, nLeaf);
[SourceClusters, SourcePerm] = sort_cluster_tree(SourceTree);
xs = xs(SourcePerm,:);
ReceiverTree = build_cluster_tree(xr, nLeaf);
[ReceiverClusters, ReceiverPerm] = sort_cluster_tree(ReceiverTree);
xr = xr(SourcePerm,:);

%% build block tree
[B_near, B_far] = build_block_tree_2(ReceiverTree, SourceTree);

%% visualise far field block tree
for b = 1 : size(B_far,1)
    i = ReceiverClusters(B_far(b,1),:);
    j = SourceClusters(B_far(b,2),:);
    patch(j(1) + j(2)*[0 0 1 1]', i(1) + i(2)*[0 1 1 0]', b*ones(4,1));
end

%% excitation
input = rand(N,1);

%% compute response with ACA
output = aca(xs, xr, int32(SourceClusters), int32(ReceiverClusters), int32(B_far-1), input);

%% compute response with direct evaluation
output_direct = zeros(N,1);
for b = 1 : length(B_far)
    i = ReceiverClusters(B_far(b,1),:);
    i = i(1) + (1:i(2));
    j = SourceClusters(B_far(b,2),:);
    j = j(1) + (1:j(2));
    
    M = 1./sqrt(...
		bsxfun(@minus, xr(i,1), xs(j,1)').^2 +...
        bsxfun(@minus, xr(i,2), xs(j,2)').^2 +...
        bsxfun(@minus, xr(i,3), xs(j,3)').^2);
    output_direct(i) = output_direct(i) + M * input(j);
end

%% compute error norm
err = log10(norm(output-output_direct) / norm(output_direct));
fprintf(1, 'Error norm: %g\n', err);
