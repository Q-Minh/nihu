%% compile C++ mex code
mex -v CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" Boonen13.mex.cpp -I../../../eigen -output Boonen13

%% Build a good little mesh
mesh = create_sphere_boundary(1,8);
k = 1;

%% Galerkin matrix C++ version
tic;
[nodes, elements] = extract_bem_mesh(mesh);
G = Boonen13(nodes, elements, k);
tHacc = toc;

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
[Hold, Gold] = bemHG(mesh, k, 'const');
tbemHG = toc;

%%
%Gold = Gold - diag(diag(Gold));
pcolor(log10(abs((G - Gold)./Gold)));
shading flat;
set(gca, 'ydir', 'reverse');
xlim([1 50]);
ylim([1 50]);
c = colorbar;
ylabel(c, 'log10 relative error');

%%
%Hold = Hold - diag(diag(Hold));
pcolor(log10(abs((G - Hold)./Hold)));
shading flat;
set(gca, 'ydir', 'reverse');
c = colorbar;
ylabel(c, 'log10 relative error');
