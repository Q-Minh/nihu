%% compile C++ mex code
mex -v CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" Boonen13.mex.cpp -I../../../eigen -output Boonen13

%% Build a good little mesh
mesh = create_sphere_boundary(1,10);
k = 1;

%% Galerkin matrix C++ version
tic;
[nodes, elements] = extract_bem_mesh(mesh);
H = Boonen13(nodes, elements, k);
tH = toc;

%% Galerkin matrix Matlab version
tic;
[x, ~, w, i] = geo2gauss(mesh, 5);
x0 = centnorm(mesh);
r = sqrt(bsxfun(@minus, x0(:,1), x(:,1)').^2 + ...
    bsxfun(@minus, x0(:,2), x(:,2)').^2 + ...
    bsxfun(@minus, x0(:,3), x(:,3)').^2);
g = exp(-1i*k*r) ./ (4*pi*r);

W = sparse(1:length(i), i, w);
GMat = g * W;
tM = toc;

%% bemHG matrix matlab version
tic;
[Hold, ~] = bemHG(mesh, k, 'const');
tbemHG = toc;

%%
pcolor(log10(abs((H - Hold)./Hold)));
shading interp;
set(gca, 'ydir', 'reverse');
c = colorbar;
ylabel(c, 'log10 relative error');
