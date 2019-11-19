%% Benchmark plane wave scattering by a rigid sphere
% Comparison with analytical results

% Relies on application:
%   helmholtz_3d_coll_mex - Conventional 3D Helmholtz BEM, collocational
%           formalism, with MEX interface

% References:
% http://ansol.us/Products/Coustyx/Validation/MultiDomain/Scattering/ ...
%   PlaneWave/HardSphere/Downloads/dataset_description.pdf

clear;

r0 = 0.4;           % Radius of scatterer sphere
n0 = 12;            % Number of elems along radius
l0 = 1;

%// scatterer mesh, a small sphere
scatterer = create_sphere_boundary(r0, n0);  

%// field point mesh, a bigger sphere
field = create_slab(l0*[-1 0 -1; 1 0 -1; 1 0 1; -1 0 1], [100 100]);

[~, keep] = mesh_select(field, sprintf('R > %.3f', r0*1.01), 'ind', 'all');
field.Elements = field.Elements(keep, :);

%// extract the mesh description matrices
[r_nodes, r_elem] = extract_core_mesh(scatterer);
[f_nodes, f_elem] = extract_core_mesh(field);

%// call C++ code at some frequency
k = min(mesh_kmax(scatterer, 9));
tic;
fprintf('Assembly at k = %.1f rad/m ... ', k);
[Ls, Ms, Lf, Mf] = helmholtz_3d_coll_mex(r_nodes, r_elem, f_nodes, f_elem, k);
fprintf('Ready in %.2f seconds.\n', toc)

%// define the incident velocity field
dir = [1 0 0];
[r_cent, r_norm] = centnorm(scatterer);
[ps_inc, qs_inc] = incident('plane', dir, r_cent, r_norm, k);

%% Compute incident field
%// analytical solution on the field
f_cent = centnorm(field);
pf_inc = incident('plane', dir, f_cent, [], k);

%// acoustic pressure on the surface and in the field points
qs_scat = -qs_inc;
fprintf('Solving BEM equation ... '); tic;
ps_scat_num = Ms \ (Ls * qs_scat);
fprintf('Ready in %.2f seconds.\n', toc);

fprintf('Computing field ... '); tic;
pf_num = Mf*ps_scat_num - Lf*qs_scat;
fprintf('Ready in %.2f seconds.\n', toc);

%% Analytical solution using spherical series expansion
order = 20;                                     % series expansion order
r_scat = sqrt(dot(r_cent, r_cent, 2));          % r on scatterer
el_scat = acos(r_cent(:,1) ./ r_scat) - pi;     % elevation on scatterer

r_fld = sqrt(dot(f_cent, f_cent, 2));           % r on field mesh
el_fld = acos(f_cent(:,1) ./ r_fld) - pi;       % elevation on field mesh

n_elem = size(r_cent, 1);                       % number of scatterer elems
n_field = size(r_fld, 1);                       % number of field elems

% Analytical series expansion of incident field on the surface
ps_inc_ana_ser = zeros(n_elem, 1);
% Analytical expression of scattered field on the surface
ps_scat_ana = zeros(n_elem, 1);
% Analytical expression of scattered field on the field mesh
pf_scat_ana = zeros(n_field, 1);

% Helper spherical functions (not included in Matlab)
% Spherical bessel function of the first kind j1
sph_j1 = @(n, z)(sqrt(pi ./ (2*z)).*besselj(n+1/2,z,1));
% Spherical hankel function of the second kind h2
sph_h2 = @(n, z)(sqrt(pi ./ (2*z)).*besselh(n+1/2,2,z));

fprintf('Calculating analytical solution ... '); tic;
for l = 0 : 1 : order
    % Precompute legendre functions
    Pl_rad = legendre(l, cos(el_scat), 'sch');
    Pl_rad = Pl_rad(1, :).';
    
    Pl_fld = legendre(l, cos(el_fld), 'sch');
    Pl_fld = Pl_fld(1, :).';
    
    % Incident field
    ps_inc_ana_ser = ps_inc_ana_ser + (2*l+1)*(1i)^(l) * sph_j1(l, k*r_scat) .* Pl_rad;

    % Scattered field coefficient
    Al = -(2*l+1)*(1i)^l * (l * sph_j1(l-1, k*r0) - (l+1) * sph_j1(l+1, k*r0)) ...
                         / (l * sph_h2(l-1, k*r0) - (l+1) * sph_h2(l+1, k*r0));
    
    % Scattered field on surface - analytical
    ps_scat_ana = ps_scat_ana + Al * sph_h2(l, k*r_scat) .* Pl_rad;
    % Scattered field on field mesh - analytical
    pf_scat_ana = pf_scat_ana + Al * sph_h2(l, k*r_fld) .* Pl_fld;
end
fprintf('Ready (in %.2f seconds)\n', toc);

%% Compute errors
err_ser = norm(ps_inc - ps_inc_ana_ser) / norm(ps_inc);
err_surf = norm(ps_scat_ana - ps_scat_num) / norm(ps_scat_ana);
err_field = norm(pf_scat_ana - pf_num) / norm(pf_scat_ana);

fprintf('Series truncation error:  %g\n', err_ser);
fprintf('Surface relative error:   %g\n', err_surf);
fprintf('Field relative error:     %g\n', err_field);

%% Plot the results on surface
fig = figure;
plot_mesh(field, abs(pf_num + pf_inc)); shading flat;
hold on;
plot_mesh(scatterer, abs(ps_scat_num + ps_inc)); view(2); axis equal;
view([0 0]);
set(gca, 'XLim', [-1 1], 'ZLim', [-1 1]);
cb = colorbar;
caxis([0 2]);
xlabel(gca, 'x [m]');
zlabel(gca, 'z [m]');
ylabel(cb, 'Sound pressure field |p| [Pa]', 'FontSize', 10);
title(sprintf('Total field at kr_0 = %.2f', k*r0));
print(fig, '-dpng', '-r100', 'bench_plane_wave_sphere.png');