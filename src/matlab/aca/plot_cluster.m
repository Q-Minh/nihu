function plot_cluster(Cluster, x, Col)

x0 = x(Cluster.ind,:);
plot3(x0(:,1), x0(:,2), x0(:,3), '.', 'Color', Col);

end