clear;

%// sphere radiator with radius R = 1, divided into 8 elements along the radius.
radiator = create_sphere_boundary(1, 30);

%// field point mesh, a line revolved around the z axis
%field = revolve_mesh(create_line([1.125, 0, 0; 3.625, 0, 0], 40), pi/100, 50, [0 0 1]);

%%
helmholtz_3d_hffmm_mex('init');

[r_nodes, r_elems] = extract_core_mesh(radiator, 'surface');
%[f_nodes, f_elems] = extract_core_mesh(field);
%%
accuracy = 3.0;
k = .8 * min(mesh_kmax(radiator));
helmholtz_3d_hffmm_mex('set', ...
    'accuracy', 3.0, 'wave_number', k);

%%
helmholtz_3d_hffmm_mex('mesh', r_nodes, r_elems);


%%
leaf_diameter = 1 / k;
helmholtz_3d_hffmm_mex('tree', 'divide_diameter', leaf_diameter);
%%
fprintf('Assembling FMM matrices ... '); tic;
helmholtz_3d_hffmm_mex('matrix');
fprintf('Ready in %.3f seconds\n', toc);

%% Prepare excitation and right hand side vector
x_src = [.2 .3 .6];
[centers, normals] = centnorm(radiator);

[p_ana, q] = incident('point', x_src, centers, normals, k);

rhs = helmholtz_3d_hffmm_mex('mvp_slp', q);

%% Prepare solver 
fprintf('Solving FMBEM system ... '); tic;
Afun = @(x)helmholtz_3d_hffmm_mex('mvp_dlp', ensure_complex(x));
[p_bem, flag, relres, iter, resvec] = gmres(Afun, rhs, [], 1e-8, 1000);
fprintf('Ready in %.3f seconds.\n', toc);

%% Evaluate and plot results
figure;
subplot(1,2,1);
plot_mesh(radiator, real(p_bem))
colorbar('Location', 'SouthOutside');
axis equal;
subplot(1,2,2);
plot_mesh(radiator, real(p_ana))
axis equal;
colorbar('Location', 'SouthOutside');
err = norm(p_ana-p_bem)/norm(p_ana);
fprintf('Relative error on surface: %.2f %%\n', err*100);

%%
helmholtz_3d_hffmm_mex('cleanup');

%%
%clear mex;