% test_telles.m
% Test Telles quadrature transform for singular and nearly singular
% integrals

% Some examples from:
% L. Jun, G. Beer, J.L. Meek: Efficient evaluation of integrals of order 
%   1/r 1/r^2, 1/r^3 using Gauss quadrature. Engineering Analysis Vol. 2,
%   Issue 3, pp. 118-128 (1985) DOI: 10.1016/0264-682X(85)90014-0

%% Test1 - singular line integral with logarithmic singularity
n_gauss = 10;
[xi, w] = gaussquad(n_gauss);

f = @(x)log(abs(.3 + x));

xi_0 = -0.3;
[xi_t, w_t] = telles_transform(xi, w, xi_0);

I_ana = -1.908598917;
I_gau = w.' * f(xi);
I_tel = w_t.' * f(xi_t);

fprintf('Test 1: int(log |0.3 + x|)dx, x = -1 .. 1\n');
fprintf('\tAnalytical: %.12f\n', I_ana);
fprintf('\tGauss  (%d points): %.12f (rel. err.: %g)\n', ...
    n_gauss, I_gau, (I_gau - I_ana) / abs(I_ana));
fprintf('\tTelles (%d points): %.12f (rel. err.: %g)\n', ...
    n_gauss, I_tel, (I_tel - I_ana) / abs(I_ana));
fprintf('\n');

%% Test2 - nearly singular line integral with 1/r^2 singularity
n_gauss = 10;
[xi, w] = gaussquad(n_gauss);

f = @(x)(1./(1.1 - x).^2);

xi_0 = 1.1;
[xi_t, w_t] = telles_transform(xi, w, xi_0);

I_ana = 9.52380952;
I_gau = w.' * f(xi);
I_tel = w_t.' * f(xi_t);

fprintf('Test 2: int(1/(1.1 - x)^2 dx, x = -1 .. 1\n');
fprintf('\tAnalytical: %.12f\n', I_ana);
fprintf('\tGauss  (%d points): %.12f (rel. err.: %g)\n', ...
    n_gauss, I_gau, (I_gau - I_ana) / abs(I_ana));
fprintf('\tTelles (%d points): %.12f (rel. err.: %g)\n', ...
    n_gauss, I_tel, (I_tel - I_ana) / abs(I_ana));
fprintf('\n');

%% Test3 - nearly singular line integral with 1/r^2 singularity

n_gauss = 10;
[xi, w] = gaussquad(n_gauss);

f = @(x)(1./(1.004 - x).^2);

xi_0 = 1.004;
[xi_t, w_t] = telles_transform(xi, w, xi_0);

I_ana = 249.500998;
I_gau = w.' * f(xi);
I_tel = w_t.' * f(xi_t);

fprintf('Test 3: int(1/(1.004 - x)^2 dx, x = -1 .. 1\n');
fprintf('\tAnalytical: %.12f\n', I_ana);
fprintf('\tGauss  (%d points): %.12f (rel. err.: %g)\n', ...
    n_gauss, I_gau, (I_gau - I_ana) / abs(I_ana));
fprintf('\tTelles (%d points): %.12f (rel. err.: %g)\n', ...
    n_gauss, I_tel, (I_tel - I_ana) / abs(I_ana));
fprintf('\n');

%% Test4 - nearly singular quadrilateral integral with 1/r singularity
n_gauss = 6;
[xi, w] = gaussquad(n_gauss);
[x, y] = ndgrid(xi, xi);
[wx, wy] = ndgrid(w, w);
xi = [x(:) y(:)];
w = wx(:) .* wy(:);
n_gauss = size(xi, 1);

f = @(x)(1./sqrt((1.004 - x(:,1)).^2 + (1.004 - x(:,2)).^2));
xi_0 = [1.004, 1.004];
[xi_t, w_t] = telles_transform(xi, w, xi_0);
I_ana = 3.476318;
I_gau = w.' * f(xi);
I_tel = w_t.' * f(xi_t);

fprintf('Test 4: int(1/sqrt((1.004-x)^2 + (1.004-y)^2) dxdy, x,y = -1 .. 1\n');
fprintf('\tAnalytical: %.12f\n', I_ana);
fprintf('\tGauss  (%d points): %.12f (rel. err.: %g)\n', ...
    n_gauss, I_gau, (I_gau - I_ana) / abs(I_ana));
fprintf('\tTelles (%d points): %.12f (rel. err.: %g)\n', ...
    n_gauss, I_tel, (I_tel - I_ana) / abs(I_ana));
fprintf('\n');

%% Test5 - 1/r in nearly singular over unit quad
n_gauss = 6;
[xi, w] = gaussquad(n_gauss);
[x, y] = ndgrid(xi, xi);
[wx, wy] = ndgrid(w, w);
xi = [x(:) y(:)];
w = wx(:) .* wy(:);
n_gauss = size(xi, 1);

f = @(x)(1./sqrt((1.02-x(:,1)).^2 + (1.02 - x(:,2)).^2));
xi_0 = [1.02, 1.02];
[xi_t, w_t] = telles_transform(xi, w, xi_0);

I_ana = 3.343673397;
I_gau = w.' * f(xi);
I_tel = w_t.' * f(xi_t);

fprintf('Test 5: int(1/sqrt((1.02-x)^2 + (1.02-y)^2) dxdy, x,y = -1 .. 1\n');
fprintf('\tAnalytical: %.12f\n', I_ana);
fprintf('\tGauss  (%d points): %.12f (rel. err.: %g)\n', ...
    n_gauss, I_gau, (I_gau - I_ana) / abs(I_ana));
fprintf('\tTelles (%d points): %.12f (rel. err.: %g)\n', ...
    n_gauss, I_tel, (I_tel - I_ana) / abs(I_ana));
fprintf('\n');

%% Test6 - 1/r^3 in nearly singular over unit quad
n_gauss = 8;
[xi, w] = gaussquad(n_gauss);
[x, y] = ndgrid(xi, xi);
[wx, wy] = ndgrid(w, w);
xi = [x(:) y(:)];
w = wx(:) .* wy(:);
n_gauss = size(xi, 1);

f = @(x)(1./((1.2-x(:,1)).^2 + (1.2 - x(:,2)).^2).^(3/2));
xi_0 = [1.2, 1.2];
[xi_t, w_t] = telles_transform(xi, w, xi_0, 0, 3);

I_ana = 2.327344790;
I_gau = w.' * f(xi);
I_tel = w_t.' * f(xi_t);

fprintf('Test 6: int(1/((1.2-x)^2 + (1.2-y)^2)^(3/2) dxdy, x,y = -1 .. 1\n');
fprintf('\tAnalytical: %.12f\n', I_ana);
fprintf('\tGauss  (%d points): %.12f (rel. err.: %g)\n', ...
    n_gauss, I_gau, (I_gau - I_ana) / abs(I_ana));
fprintf('\tTelles (%d points): %.12f (rel. err.: %g)\n', ...
    n_gauss, I_tel, (I_tel - I_ana) / abs(I_ana));
fprintf('\n');

%% Test7 - 1/r^3 in nearly singular at some distance
n_gauss = 8;
[xi, w] = gaussquad(n_gauss);
[x, y] = ndgrid(xi, xi);
[wx, wy] = ndgrid(w, w);
xi = [x(:) y(:)];
w = wx(:) .* wy(:);
n_gauss = size(xi, 1);

zeta = .25;
f = @(x)(1./((0.7-x(:,1)).^2 + (0.5 - x(:,2)).^2 + zeta^2).^(3/2));
xi_0 = [0.7, 0.5];
[xi_t, w_t] = telles_transform(xi, w, xi_0, zeta, 3);

[XI, W] = gaussquad2(100);
I_ana = W.' * f(XI);
I_gau = w.' * f(xi);
I_tel = w_t.' * f(xi_t);

fprintf('Test 7: int(1/((0.7-x)^2 + (0.5-y)^2 + 0.25^2)^(3/2) dxdy, x,y = -1 .. 1\n');
fprintf('\tBrute force quadrature: %.12f\n', I_ana);
fprintf('\tGauss  (%d points): %.12f (rel. err.: %g)\n', ...
    n_gauss, I_gau, (I_gau - I_ana) / abs(I_ana));
fprintf('\tTelles (%d points): %.12f (rel. err.: %g)\n', ...
    n_gauss, I_tel, (I_tel - I_ana) / abs(I_ana));
fprintf('\n');

%% Test8 - laplace_3d_HSP_kernel in nearly singular at some distance
clear;

n_gauss = 18;
[xi, w] = gaussquad2(2*n_gauss-1, 4);
n_gauss = size(xi, 1);

coords = [
    0 0 0
    1 0 0
    1 .5 0
    0 1 .2
    ];
x = [.5 .7 .2];

nx = [0 0 1];

lset = ShapeSet.LinearQuad;
nset = ShapeSet.LinearQuad;

xi0 = be_invmap(coords.', x, lset);
zeta = xi0(3);
xi0 = xi0(1:2);

[xi_t, w_t] = telles_transform(xi, w, xi0, zeta);

figure;
plot3(xi(:,1), xi(:,2), w, '.', ...
    xi_t(:,1), xi_t(:,2), w_t, '.');

[L, dL] = lset.eval(xi_t);
y = L * coords;
jvec = cross(dL(:,:,1) * coords, dL(:,:,2) * coords, 2);
jac = sqrt(dot(jvec, jvec, 2));
ny = bsxfun(@times, jvec, 1./jac);

addpath(fullfile(pwd(), '../helmholtz_bem_implementation'));

f = laplace_3d_hsp_kernel(x, nx, y, ny);
N = nset.eval(xi_t);

I_tel = (w_t.*jac).' * bsxfun(@times, f, N);

fprintf('\tTelles (%d points): \n%.12f\n%.12f\n%.12f\n%.12f\n', ...
    n_gauss, I_tel);
%% Test 9 - Tria
order = 1000;
n_gauss = ceil((order+1)/2);
[xi, w] = duffy_gaussian_quadrature(Domain.Tria, order, 1);

n_telles = 20;
xi_0 = [-0.15, 0.15];
[xi_t, w_t] = telles_tria(n_telles, xi_0, 0);


f = @(x)(1./sqrt((xi_0(1)-x(:,1)).^3 + (xi_0(2) - x(:,2)).^2));

I_duffy = w.' * f(xi);
I_telles = w_t.' * f(xi_t);

fprintf('Test 9: 1/r^2 over a triangle\n');
fprintf('\tDuffy  (%d points):\t%.12f\n', size(xi, 1), I_duffy);
fprintf('\tTelles (%d points):\t%.12f (rel.diff.: %.3g)\n', size(xi_t, 1), I_telles, (I_telles-I_duffy)/abs(I_duffy));

