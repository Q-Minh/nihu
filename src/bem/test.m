make

%%
mesh = create_sphere_boundary(1,20);
field = create_slab([2 0 0; 4 2 0], [50, 50]);

k = mesh_kmax(mesh, 7) * .8;

[nodes, elements] = extract_bem_mesh(mesh);
points = field.Nodes(:,2:4);

tic;
[GB, HB] = Boonen13(nodes, elements, k, points);
tBoonen_acc = toc;

tic;
[H, G] = bemHG(mesh, k, 'lin', points);
tOldSchool = toc;

%%
figure;
subplot(1,2,1);
pcolor(log10(abs(deltaG)));
shading interp
subplot(1,2,2);
pcolor(log10(abs(deltaH)));
shading interp

%%
v = mesh.Nodes(:,2);
subplot(1,3,1);
plot_mesh(field, real(G*v));
shading interp;
C = caxis;
subplot(1,3,2);
plot_mesh(field, real(B*v));
shading interp;
caxis(C);
subplot(1,3,3);
plot_mesh(field, log10(abs((G*v-B*v)./(G*v))));
shading interp;
colorbar;