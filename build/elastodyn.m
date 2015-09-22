Us = import_matrix('Us.mtx');
Uf = import_matrix('Uf.mtx');
Ts = import_matrix('Ts.mtx');
Tf = import_matrix('Tf.mtx');

mesh0 = import_mesh('../src/test/stallone_unit/sphere_1_5.off');
mesh = import_mesh('../src/test/stallone_unit/sphere_2_5.off');
xs = centnorm(mesh0);
xf = centnorm(mesh);

R0 = mean(sqrt(dot(xs,xs,2)));
R = mean(sqrt(dot(xf,xf,2)));

t = reshape(xs', numel(xs),1);
u0 = (Ts - .5*eye(size(Ts))) \ (Us * t);
u = Tf * u0 - Uf * t;
u0vec = reshape(u0, 3, []).';
uvec = reshape(u, 3, []).';
u0a = sqrt(dot(u0vec, u0vec, 2));
ua = sqrt(dot(uvec, uvec, 2));

nu = .33; rho = 100; mu = 1e8;
freq = 50;
g = elastodynamic_pulsating_sphere(R0, R, mu, rho, nu, freq);
g0 = elastodynamic_pulsating_sphere(R0, R0, mu, rho, nu, freq);
