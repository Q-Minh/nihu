[pot_bb, pot_path] = read_epspath('potato.eps');
siz = norm(diff(pot_bb));
Le = siz/100;
radiator = meshpath(pot_path, Le);
radiator = translate_mesh(scale_mesh(radiator, 1/siz), [.2 .2]);

field = create_slab([1.5 1], [150 100]);

[r_nodes, r_elem] = extract_core_mesh(radiator);
[f_nodes, f_elem] = extract_core_mesh(field);
f_elem(f_elem(:,1) == 32404,1) = 22404;

k = .7 * min(min(mesh_kmax(radiator)), min(mesh_kmax(field)));
[Ls, Ms, Lf, Mf] = helmholtz_bem_2d(r_nodes, r_elem, f_nodes, f_elem, k);

x0 = [0. 0. 0.];
r_cent = centnorm(radiator);
ps = incident('line', x0, r_cent, [], k);

f_cent = centnorm(field);
pf_ana = incident('line', x0, f_cent, [], k);

I = speye(size(Ms));
qs_int = Ls \ (Ms - .5 * I) * ps;
qs_ext = Ls \ (Ms + .5 * I) * ps;

pf = Lf * (qs_ext - qs_int);

figure;
hold on;
plot(x0(:,1), x0(:,2), 'k*');
plot_mesh(field, real(pf));
shading flat;
plot_mesh(radiator);
axis equal tight;
c1 = colorbar;
ylabel(c1, 'Real part of numerical solution');
