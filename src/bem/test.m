make

%%
mesh = create_sphere_boundary(1,20);
% mesh = quad2tria(mesh);
field = create_slab([2 0 0; 4 2 0], [20, 20]);

k = min(mesh_kmax(mesh, 7)) * .5;

[nodes, elements] = extract_bem_mesh(mesh);
points = field.Nodes(:,2:4);

%% matrices
tic;
[GB, HB] = Boonen13(nodes, elements, k, points);
tBoonen = toc;

tic;
[H0, G0] = bemHG(mesh, k, 'lin', points);
tOldSchool = toc;

%% excitation and response
r0 = [0 0 0];
[ps, qs] = incident('point', r0, nodes, nodes, k);

pf0 = incident('point', r0, points, [], k);
pf_Boonen = HB * ps - GB * qs;
pf_OldSchool = H0 * ps - G0 * qs;

%%
figure(1);

subplot(2,2,1);
plot_mesh(field, pf0);
shading interp

subplot(2,2,2);
plot_mesh(field, pf_Boonen);
shading interp

subplot(2,2,3);
plot_mesh(field, log10(abs(pf_Boonen./pf0-1)));
colorbar;
shading interp

subplot(2,2,4);
plot_mesh(field, log10(abs(pf_OldSchool./pf0-1)));
colorbar;
shading interp
