%% create surface and field point meshes
Ne = 5;
mesh = create_sphere_boundary(1, Ne);
field = create_sphere_boundary(2, 5);

mater.mu = 1;
mater.nu = .33;

% surface element centers and field points
[xs, ns] = centnorm(mesh);
xf = centnorm(field);

%% compute BEM system matrices
tic;
% preallocate
M = zeros(3*size(xf,1), 3*size(xs,1));
L = zeros(3*size(xf,1), 3*size(xs,1));
% fill row-by-row
for n = 1 : size(xf,1)
    progbar(1, size(xf, 1), n);
    x0 = xf(n,1:3); % collocation point
    rowind = (n-1)*3+(1:3);
    for e = 1 : size(xs, 1)
        colind = (e-1)*3+(1:3);
        coords = mesh.Nodes(mesh.Elements(e, 5:7), 2:4);
        L(rowind,colind) = elasto_3d_regular('U', mater, x0, coords);
        M(rowind,colind) = elasto_3d_regular('T', mater, x0, coords);
    end
end
toc;

figure;
surf(L(1:3:end, 1:3:end));
figure;
surf(M(1:3:end, 1:3:end));

%% BC
% source point
x0 = [.2 .2 0];

% potential and normal derivative on surface
us = zeros(size(L,2), 1);
ts = zeros(size(L,2), 1);
for e = 1 : size(xs,1)
    us(3*(e-1)+(1:3)) = Ukernel(mater, x0, xs(e,:)) * [0; 0; 1];
    ts(3*(e-1)+(1:3)) = Tkernel(mater, x0, xs(e,:), ns(e,:)) * [0; 0; 1];
end

% potential and on field
uf_anal = zeros(size(L,1), 1);
for f = 1 : size(xf,1)
    uf_anal(3*(f-1)+(1:3)) = Ukernel(mater, x0, xf(f,:)) * [0; 0; 1];
end

uf = M*us - L*ts;

result = sqrt(dot(reshape(uf, 3, []), reshape(uf, 3, []), 1)).';
result_anal = sqrt(dot(reshape(uf_anal, 3, []), reshape(uf_anal, 3, []), 1)).';

figure;
plot_mesh(field, result);
caxis([2 5]*1e-10);
figure;
plot_mesh(field, result_anal);
caxis([2 5]*1e-10);

error = log10(abs(result./result_anal-1));

fprintf(1, 'Minimal log10 error: %g\n', min(error));
fprintf(1, 'Mean log10 error: %g\n', mean(error));
fprintf(1, 'Maximal log10 error: %g\n', max(error));
