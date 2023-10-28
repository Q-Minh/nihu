clear;

%%
dim = 2;

msh = create_circle(1, 70);
xc = centnorm(msh);
xc = xc(:,1:dim);

%% generate tree with normal near field
T = build_regular_clustertree(6, xc);

figure;
plot_cluster(T(6), 192);
set(gcf, 'renderer', 'painters');

%%
nearfun = @(d)(d(:,end)-1).^2 - dot(d(:,1:end-1), d(:,1:end-1), 2) < 4-1e-2;

T = build_regular_clustertree(6, xc, [], nearfun);

figure;
plot_cluster(T(6), 192);
set(gcf, 'renderer', 'painters');


%%
