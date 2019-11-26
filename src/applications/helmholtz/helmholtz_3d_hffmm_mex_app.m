clear;

mex_fun = @helmholtz_3d_hffmm_mex;

%// sphere radiator with radius R = 1, divided into 30 elements along the radius.
radiator = create_sphere_boundary(1, 30);

%// field point mesh, a line revolved around the z axis
field = revolve_mesh(create_line([1.125, 0, 0; 5.625, 0, 0], 180), pi/200, 100, [0 0 1]);

%% Initialize the FMM and set parameters
mex_fun('init');
[r_nodes, r_elems] = extract_core_mesh(radiator, 'surface');

k = .8 * min(mesh_kmax(radiator));
mex_fun('set', ...
    'accuracy', 3.0, ...
    'far_field_order', 5, ...
    'wave_number', k);

%% Assign the mesh and create the cluster tree
mex_fun('mesh', r_nodes, r_elems);
mex_fun('tree', 'divide_diameter', 1/k);
%% Assemble the FMM matrices for the surface
fprintf('Assembling FMM matrices ...\n'); tic;
mex_fun('matrix');
fprintf('Ready in %.3f seconds\n', toc);

%% Prepare excitation and right hand side vector
x_src = [.2 .3 .6];
[centers, normals] = centnorm(radiator);
[ps_ana, qs] = incident('point', x_src, centers, normals, k);
rhs = mex_fun('mvp_slp', qs);

%% Prepare solver 
fprintf('Solving FMBEM system ... '); tic;
Afun = @(x)mex_fun('mvp_dlp', ensure_complex(x));
[ps_bem, flag, relres, iter, resvec] = gmres(Afun, rhs, [], 1e-8, 1000);
fprintf('Ready in %.3f seconds.\n', toc);

%% Evaluate and plot surface results
figure;
subplot(1,2,1);
plot_mesh(radiator, real(ps_bem))
colorbar('Location', 'SouthOutside');
axis equal;
subplot(1,2,2);
plot_mesh(radiator, real(ps_ana))
axis equal;
colorbar('Location', 'SouthOutside');
err = norm(ps_ana-ps_bem)/norm(ps_ana);
fprintf('Relative error on surface: %.2f %%\n', err*100);

%% Plot iteration results
figure;
plot(log10(resvec));
xlabel('#iterations');
ylabel('log10 relative residual');

%% Clean up the surface computations
mex_fun('cleanup');

%% Initialize for field calculations
mex_fun('init');
[f_nodes, f_elems] = extract_core_mesh(field, 'surface');

mex_fun('set', ...
    'accuracy', 3.0, ...
    'far_field_order', 5, ...
    'wave_number', k);

%% Assign the mesh and create the cluster tree
mex_fun('mesh', r_nodes, r_elems, f_nodes, f_elems);
mex_fun('tree', 'divide_diameter', 1/k);
%% Assemble the FMM matrices for the surface
fprintf('Assembling FMM matrices ...\n'); tic;
mex_fun('matrix');
fprintf('Ready in %.3f seconds\n', toc);

%% Compute analytical solution
centers = centnorm(field);
pf_ana = incident('point', x_src, centers, [], k);

%% Compute FMBEM solution
pf_bem = mex_fun('mvp_dlp', ensure_complex(ps_bem)) ...
       - mex_fun('mvp_slp', ensure_complex(qs));
   
%%
figure;
subplot(1,2,1);
plot_mesh(radiator, real(ps_bem));
hold on
plot_mesh(field, real(pf_bem));
colorbar('Location', 'SouthOutside');
axis equal;
subplot(1,2,2);
plot_mesh(radiator, real(ps_ana));
hold on
plot_mesh(field, real(pf_ana));
axis equal;
colorbar('Location', 'SouthOutside');
%% Evaluate error
err = norm(pf_ana-pf_bem)/norm(pf_ana);
fprintf('Relative error on field: %.2f %%\n', err*100);
%% Clean up
mex_fun('cleanup');