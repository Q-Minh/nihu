mex -v CXXFLAGS="\$CXXFLAGS -std=c++0x -O3" Boonen13.mex.cpp -I../../../eigen -output Boonen13


%% Build a good little mesh
mesh = create_sphere_boundary(1,5);
k = 1;

%% Galerkin matrix C++ version
tic;
[nodes, elements] = extract_bem_mesh(mesh);
G = Boonen13(nodes, elements, k);
tG = toc;

%% Galerkin matrix Matlab version
tic;
[x, ~, w, i] = geo2gauss(mesh, 5);
r = sqrt(bsxfun(@minus, x(:,1), x(:,1)').^2 + ...
    bsxfun(@minus, x(:,2), x(:,2)').^2 + ...
    bsxfun(@minus, x(:,3), x(:,3)').^2);
g = exp(-1i*k*r) ./ (4*pi*r);

W = sparse(1:length(i), i, w);
GMat = W' * g * W;
tM = toc;

%%
pcolor(log10(abs((G - GMat)./GMat)));
shading interp;
set(gca, 'ydir', 'reverse');
c = colorbar;
ylabel(c, 'log10 relative error');
