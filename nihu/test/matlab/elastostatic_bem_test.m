clear
close all;

%% create meshes
Rs = 1;
Rf = 3;
N = 6;

radiator = create_sphere_boundary(Rs, N);
field = create_slab([10 .1], [100, 1]);
field = translate_mesh(field, [1.1*Rs-5e-2 -5e-2 0]);

[xs, ns] = centnorm(radiator);
xf = centnorm(field);

%% material properties and frequency
nu = .3;    % [-]
mu = 1e8;   % Pa

%% compute system matrices
[r_nodes, r_elements] = extract_core_mesh(radiator);
[f_nodes, f_elements] = extract_core_mesh(field);

tic;
[Us, Ts, Uf, Tf] = elastostatic_bem(r_nodes, r_elements, f_nodes, f_elements,...
    nu, mu);
t_mat = toc;


%% Pulsating test
tsvec = ns;
ts = reshape(tsvec.', [], 1);
us = (Ts - .5*eye(size(Ts))) \ (Us * ts);
uf = Tf * us - Uf * ts;

usvec = reshape(us, 3, []).';
ufvec = reshape(uf, 3, []).';

usn = dot(usvec, ns, 2);
ufn = ufvec(:,1);

R = xf(:,1);
gs = elastodynamic_pulsating_sphere(Rs, Rs, mu, 1, nu, 0);
gf = elastodynamic_pulsating_sphere(Rs, R, mu, 1, nu, 0);

figure;
plot(xf(:,1), [ufn gf]);
legend('numerical', 'analytical');
hold on;
plot(Rs, [usn(1), gs], '*');

%% Transparent test
load_normal = [1 0 0].';
usvec = zeros(size(xs));
ufvec_anal = zeros(size(xf));
x0 = [.1 .2 -.3];
for i = 1 : size(xs,1)
    rvec = xs(i,:)-x0;
    usvec(i,:) = elastodynamic_green(rvec(:), mu, 1, nu, 0) * load_normal;
end
for i = 1 : size(xf,1)
    rvec = xf(i,:)-x0;
    ufvec_anal(i,:) = elastodynamic_green(rvec(:), mu, 1, nu, 0) * load_normal;
end

us = reshape(usvec.', [], 1);
ts = Us \ ((Ts - .5*eye(size(Ts))) * us);
uf = Tf * us - Uf * ts;
ufvec = reshape(uf, 3, []).';

d = 3;
figure;
plot(xf(:,1), [ufvec(:,d) ufvec_anal(:,d)]);
