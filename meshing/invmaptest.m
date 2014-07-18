clear;
clc;

%%
mesh = create_brick([1 1 1], [10 10 10]);
mesh.Nodes(:,4) = mesh.Nodes(:,4) + 0.6*mesh.Nodes(:,2).^2;
mesh.Nodes(:,3) = mesh.Nodes(:,3) - 0.6*mesh.Nodes(:,2).^2;
mesh.Nodes(:,2) = mesh.Nodes(:,2) + 0.6*mesh.Nodes(:,3).^2;
p = sum(mesh.Nodes(:,2:4),2).^2;

%%
plane = create_slab([2 2], [200 200]);
plane = translate_mesh(plane, mean(mesh.Nodes(:,2:4),1)-[1 1 0]);

%%
tic;
[Ap, Av] = fe_interp(mesh, plane.Nodes(:,2:4));
toc;

%%
figure;
plot_mesh(mesh, p);
alpha .2;
plot_mesh(plane, Ap*p);
shading interp;