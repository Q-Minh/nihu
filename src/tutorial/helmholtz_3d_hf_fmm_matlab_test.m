clear;

%// sphere radiator with radius R = 1, divided into 8 elements along the radius.
radiator = create_sphere_boundary(1, 10);

%// field point mesh, a line revolved around the z axis
%field = revolve_mesh(create_line([1.125, 0, 0; 3.625, 0, 0], 40), pi/100, 50, [0 0 1]);

%%
helmholtz_3d_hf_fmm_matlab('init');

[r_nodes, r_elems] = extract_core_mesh(radiator, 'surface');
%[f_nodes, f_elems] = extract_core_mesh(field);
%%
helmholtz_3d_hf_fmm_matlab('set', ...
    'accuracy', 3.0, 'wave_number', 6);

%%
helmholtz_3d_hf_fmm_matlab('mesh', r_nodes, r_elems);


%%
n_levels = 6;
helmholtz_3d_hf_fmm_matlab('tree', 'divide_depth', n_levels);
%%
fprintf('Assembling FMM matrices ... '); tic;
helmholtz_3d_hf_fmm_matlab('matrix');
fprintf('Ready in %.3f seconds\n', toc);

%%
n_elems = size(radiator.Elements, 1);
q = complex(ones(n_elems, 1));
%%
p = helmholtz_3d_hf_fmm_matlab('mvp_slp', q);

%%
helmholtz_3d_hf_fmm_matlab('cleanup');

%%
%clear mex;