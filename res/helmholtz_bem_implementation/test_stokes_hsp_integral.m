clear;

%%
hsp_kernel = @(x, nx, y, ny)laplace_3d_hsp_kernel(x, nx, y, ny);
slp_kernel = @(x,y)laplace_3d_slp_kernel(x, [], y, []);
grad_kernel = @(x,y)laplace_3d_grad_kernel(x, [], y, []);

% elem lset and nset
lset = ShapeSet.LinearQuad;
nset = ShapeSet.LinearQuad;

% elem corner nodes
corners = [
    0 0 0
    1 0 0
    1 .5 0
    0 1 .2
    ];

% singular point in physical domain
x = [.5 .5 .2];
nx = [1 0 0];

%% Local Taylor series expansion of the shape function
% determine xi0 and shape functions in xi0
xi0 = be_invmap(corners.', x, lset);
xi0 = xi0(1:2);
[N0, dN0] = nset.eval(xi0);
N0 = N0.';
dN0 = [
    dN0(:,:,1)
    dN0(:,:,2)
    ];
dN0(3,:) = 0;

[L0, dL0] = lset.eval(xi0);
y0 = L0 * corners;
dy0 = [
    dL0(:,:,1) * corners
    dL0(:,:,2) * corners
    ];
dy0(3,:) = cross(dy0(1,:), dy0(2,:));
gradN0 = dN0.' / dy0.';

%% compute with stokes
Istokes = stokes_hsp_integral(corners, lset, x, nx, ...
    y0, N0, gradN0, ...
    slp_kernel, grad_kernel);

%% compute linear term numerically
order = 60;
ncorners = size(corners,1);
[xi, w] = gaussquad2(order, ncorners);
[L, dL] = lset.eval(xi);
y = L * corners;
jvec = cross(dL(:,:,1) * corners, dL(:,:,2) * corners, 2);
jac = sqrt(dot(jvec, jvec, 2));
ny = bsxfun(@times, jvec, 1./jac);
G = hsp_kernel(x, nx, y, ny);
N = nset.eval(xi);
Nlin = bsxfun(@plus, N0, gradN0 * bsxfun(@minus, y, y0).').';
Isurf = Nlin.' * diag(w) * (G .* jac);

%%
fprintf(1, 'Error: %g\n', norm(Isurf-Istokes) / norm(Isurf));

figure;
plot3(xi(:,1), xi(:,2), bsxfun(@times, Nlin, (G .* jac)), '.');
title('Surface integral in numerical method');