clear

Rs = 1;
Rf = 3;
N = 6;

radiator = create_sphere_boundary(Rs, N);
field = create_slab([10 .1], [100, 1]);
field = translate_mesh(field, [2*Rs-5e-2 -5e-2 0]);

[r_nodes, r_elements] = extract_core_mesh(radiator);
[f_nodes, f_elements] = extract_core_mesh(field);

nu = .3;
mu = 1e8;
rho = 100;
freq = 500;

tic;
[Us, Ts, Uf, Tf] = elastodynamic_bem(r_nodes, r_elements, f_nodes, f_elements,...
    nu, rho, mu, freq);
t_mat = toc;

[xs, ns] = centnorm(radiator);
xf = centnorm(field);

tsvec = ns;
ts = reshape(tsvec.', numel(tsvec), 1);
us = (Ts - .5*eye(size(Ts))) \ (Us * ts);
uf = Tf * us - Uf * ts;

usvec = reshape(us, 3, []).';
ufvec = reshape(uf, 3, []).';

usn = dot(usvec, ns, 2);
ufn = ufvec(:,1);

R = xf(:,1);
g = elastodynamic_pulsating_sphere(Rs, R, mu, rho, nu, freq);
