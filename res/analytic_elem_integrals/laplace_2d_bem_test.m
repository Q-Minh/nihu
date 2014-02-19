clear all;

%% create surface and field point meshes
Ne =19;
mesh = create_circle_boundary(1, Ne);
field = create_circle_boundary(2, Ne);

% surface element centers and field points
[xs, ns] = centnorm(mesh);
[xf, nf] = centnorm(field);

%% compute BEM system matrices
tic;
% preallocate
M = zeros(size(xf,1), size(xs,1));
L = zeros(size(xf,1), size(xs,1));
% fill row-by-row
for n = 1 : size(xf,1)
    x0 = xf(n,1:2); % collocation point
    for e = 1 : size(xs, 1)
        % elem nodal coordinates
        y1 = mesh.Nodes(mesh.Elements(e,5),2:3);
        y2 = mesh.Nodes(mesh.Elements(e,6),2:3);
        % transform to local coordinates
        t = y1 - y2;
        t = t / norm(t);
        nor = t * [0 1; -1 0];
        xi1 = (y1 - x0) * t.';
        xi2 = (y2 - x0) * t.';
        eta = (x0-y1) * nor.';
        
        % analytical elements
        M(n,e) = laplace_2d_regular('dlp', xi1, xi2, eta, 1, 0, 0);
        L(n,e) = laplace_2d_regular('slp', xi1, xi2, eta, 1, 0, 0);
    end
end
toc;

%% BC
% source point
x0 = [.2 .2 0];

% distance to surface
rvec = bsxfun(@minus, xs, x0);
r = sqrt(dot(rvec, rvec, 2));
% potential and normal derivative on surface
us = -log(r)/2/pi;
qs = -1./r.^2 .* dot(rvec, ns, 2) / 2/pi;

% distance to field points
rvec = bsxfun(@minus, xf, x0);
r = sqrt(dot(rvec, rvec, 2));
% potential at field points
uf_anal = -log(r)/2/pi;
% BEM results
uf = M*us - L*qs;

error = log10(abs(uf./uf_anal-1));

fprintf(1, 'Maximal log10 error: %g\n', max(error));
