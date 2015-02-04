clear
close all;

Rs = 1;
Rf = 3;
N = 6;

radiator = create_sphere_boundary(Rs, N);
field = create_slab([10 .1], [100, 1]);
field = translate_mesh(field, [1.1*Rs-5e-2 -5e-2 0]);

[r_nodes, r_elements] = extract_core_mesh(radiator);
[f_nodes, f_elements] = extract_core_mesh(field);

nu = .3;
mu = 1e8;
rho = 100;
freq = 500;

tic;
[Us, Ts, Uf, Tf] = elastostatic_bem(r_nodes, r_elements, f_nodes, f_elements,...
    nu, mu);
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
gs = elastodynamic_pulsating_sphere(Rs, Rs, mu, rho, nu, 0);
gf = elastodynamic_pulsating_sphere(Rs, R, mu, rho, nu, 0);

figure;
subplot(2,1,1);
plot(xf(:,1), [real(ufn), imag(ufn), abs(ufn)]);
hold on;
plot(Rs, real(usn(1)), 'b*', Rs, imag(usn(1)), 'g*', Rs, abs(usn(1)), 'r*');
ylim([-3 3]*1e-9);
subplot(2,1,2);
plot(xf(:,1), [real(gf), imag(gf), abs(gf)]);
ylim([-3 3]*1e-9);
hold on;
plot(Rs, real(gs(1)), 'b*', Rs, imag(gs(1)), 'g*', Rs, abs(gs(1)), 'r*');

