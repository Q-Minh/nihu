clear;

%% create surface and field point meshes
Ne = 5;
mesh = create_sphere_boundary(1, Ne);
field = create_sphere_boundary(2, Ne);

mater.mu = 1;
mater.nu = .33;

% surface element centers and field points
[xs, ns] = centnorm(mesh);
xf = centnorm(field);

%% compute BEM system matrices
tic;
% preallocate
T = zeros(3*size(xf,1), 3*size(xs,1));
U = zeros(3*size(xf,1), 3*size(xs,1));
% fill row-by-row
for n = 1 : size(xf,1)
    progbar(1, size(xf, 1), n);
    x0 = xf(n,1:3); % collocation point
    rowind = (n-1)*3+(1:3);
    for e = 1 : size(xs, 1)
        colind = (e-1)*3+(1:3);
        coords = mesh.Nodes(mesh.Elements(e, 5:7), 2:4);
        U(rowind,colind) = elasto_3d_regular('U', mater, x0, coords);
        T(rowind,colind) = elasto_3d_regular('T', mater, x0, coords);
    end
end
toc;

figure;
surf(U(1:3:end, 1:3:end));
figure;
surf(T(1:3:end, 1:3:end));

%% Load Boonen system matrices
Us = importdata('Us.mtx');
Uf = importdata('Uf.mtx');
Ts = importdata('Ts.mtx');
Tf = importdata('Tf.mtx');

%% BC
% source point
x0 = [0 0 0];

% potential and normal derivative on surface
us = zeros(size(mesh.Elements,1)*3, 1);
ts = us;
for e = 1 : size(xs,1)
    us(3*(e-1)+(1:3)) = Ukernel(mater, x0, xs(e,:)) * [0; 0; 1];
    ts(3*(e-1)+(1:3)) = Tkernel(mater, x0, xs(e,:), ns(e,:)) * [0; 0; 1];
end

%% Solution with Boonen matrices
us_Boonen = (Ts - .5*eye(size(Ts))) \ (Us * ts);
anal = reshape(us,3,[]);
diff = reshape(us_Boonen,3,[]) - anal;
err_s = .5 * log10( dot(diff,diff,1)./dot(anal, anal, 1) );
figure;
dir = 2;
a1 = subplot(2,2,1);
plot_mesh(mesh, err_s);
a2 = subplot(2,2,2);
plot_mesh(mesh, us_Boonen(dir:3:end));
a3 = subplot(2,2,4);
plot_mesh(mesh, us(dir:3:end));
linkprop([a2,a3], 'CLim');
linkprop([a1,a2,a3], 'View');

%%
% potential and on field
uf_anal = zeros(size(U,1), 1);
for f = 1 : size(xf,1)
    uf_anal(3*(f-1)+(1:3)) = Ukernel(mater, x0, xf(f,:)) * [0; 0; 1];
end

uf = T*us - U*ts;

result = sqrt(dot(reshape(uf, 3, []), reshape(uf, 3, []), 1)).';
result_anal = sqrt(dot(reshape(uf_anal, 3, []), reshape(uf_anal, 3, []), 1)).';

figure;
plot_mesh(field, result);
figure;
plot_mesh(field, result_anal);

error = log10(abs(result./result_anal-1));

fprintf(1, 'Minimal log10 error: %g\n', min(error));
fprintf(1, 'Mean log10 error: %g\n', mean(error));
fprintf(1, 'Maximal log10 error: %g\n', max(error));

save results
