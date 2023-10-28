clear;

%%
R = .5;
H = 2;
Le = .1;
mesh = create_slab(2*R, ceil(2*R/Le));
Nz = ceil(H / Le);
dz = H / Nz;
mesh = extrude_mesh(mesh, [0 0 dz], Nz);
mesh = get_boundary(mesh);
mesh = flip_mesh(mesh);
mesh = drop_mesh_IDs(drop_unused_nodes(mesh));

%% clamped mesh
[~, clamp_ind] = mesh_select(mesh, 'abs(z) < 1e-4', 'ind', 'all');
free_ind = setdiff(1 : size(mesh.Elements,1), clamp_ind);
cmesh = mesh;
cmesh.Elements = cmesh.Elements(free_ind,:);
cmesh = drop_mesh_IDs(drop_unused_nodes(cmesh));

%% excitation
[~, top_ind] = mesh_select(cmesh, sprintf('abs(z - %g) < 1e-4', H), 'ind', 'all');
t_exc = zeros(size(cmesh.Elements,1), 3);
t_exc(top_ind, 1) = +1;

%%
export_off_mesh(cmesh, 'data/cylinder_10cm.off');
fid = fopen('data/pull_10cm.xct', 'wt');
fprintf(fid, '%g\n', reshape(t_exc', [], 1));
fclose(fid);

%% call NiHu here

%%
fid = fopen('data/pull_10cm.res', 'rt');
cu_res = fscanf(fid, '%g', inf);
cu_res = reshape(cu_res(:), 3, [])';
fclose(fid);

u_res = zeros(size(mesh.Elements, 1), 3);
u_res(free_ind, :) = cu_res;

[~, A] = adjacency(mesh);
A = A';
u_res_node = bsxfun(@times, A * u_res, 1./sum(A, 2));

pmesh = mesh;
mag = 1e7;
pmesh.Nodes(:,2:4) = pmesh.Nodes(:,2:4) + mag * u_res_node;

figure;
plot_mesh(pmesh)
axis equal;
