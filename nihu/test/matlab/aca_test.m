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
[B_near, B_far] = build_block_tree_2(ReceiverTree, SourceTree, .8);

%% visualise far field block tree
display_block_structure(ReceiverClusters, SourceClusters, B_far);

%% excitation
input = ones(N,1);

%% compute response with ACA
output = aca(xs, xr, int32(SourceClusters), int32(ReceiverClusters), int32(B_far-1), input);

%% compute response with direct evaluation
M = @(i,j)laplace_matrix(i,j,xr,xs);
output_direct = zeros(N,1);
for b = 1 : length(B_far)
    i = ReceiverClusters(B_far(b,1),:);
    i = i(1) + (1:i(2));
    j = SourceClusters(B_far(b,2),:);
    j = j(1) + (1:j(2));
    
    output_direct(i) = output_direct(i) + M(i,j) * input(j);
end

%% compute error norm
err = log10(norm(output-output_direct) / norm(output_direct));
fprintf(1, 'Error norm: %g\n', err);
