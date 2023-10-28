function elasto_3d_bem_test

%% create surface and field point meshes
Ne = 4;
mesh = quad2tria(create_sphere_boundary(1, Ne));
field = create_sphere_boundary(2, 3);

mater.mu = 1e8;
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

figure;
plot_mesh(field, uf(1:3:end));

figure;
plot_mesh(field, uf_anal(1:3:end));

error = log10(abs(uf./uf_anal-1));

fprintf(1, 'Maximal log10 error: %g\n', max(error));

end


function K = elasto_3d_regular(type, mater, x0, coords)
switch type
    case 'T'
        kernel = @Tkernel;
    case 'U'
        kernel = @Ukernel;
end
n = cross(coords(2,:) - coords(1,:), coords(3,:) - coords(1,:));
jac = norm(n);
n = n / jac;

[xi, w] = gaussquad2(5, 3);
xg = shapefun(xi, 23) * coords;

K = zeros(3,3);
for g = 1 : size(xg,1)
    K = K + w(g) * kernel(mater, x0, xg(g,:), n);
end

K = jac * K;
end


function U = Ukernel(mater, x, y, ~)
rvec = y(:) - x(:);
r = norm(rvec);
rdy = rvec / r;
A = 1/(16*pi*mater.mu*(1-mater.nu)*r);
U = A * ((3-4*mater.nu)*eye(3) + (rdy*rdy.'));
end

function T = Tkernel(mater, x, y, n)
rvec = y(:) - x(:);
r = norm(rvec);
rdy = rvec / r;
rdn = dot(rdy, n(:), 1);
B = 1 / (8*pi*(1-mater.nu)*r^2);
T = -B * rdn * ((1-2*mater.nu)*eye(3) + 3*(rdy*rdy.'));
for i = 1 : 3
    for k = 1 : 3
        T(i,k) = T(i,k) + (1-2*mater.nu)*B * (rdy(k)*n(i) - rdy(i)*n(k));
    end
end
end
