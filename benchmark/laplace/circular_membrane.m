% circular_membrane.m
% Circular membrane under static load at a single point

% Comparison of boundary element method and Green's function obtained using
% the method of mirror images for the circular membrane

r0 = .3;
n0 = 40;
membrane = create_circle(r0, n0);
boundary = drop_mesh_IDs(drop_unused_nodes(get_boundary(membrane)));
% Interior problem, normals are flipped
boundary = flip_mesh(boundary);

%% Extract the meshes, assemble matrices
[b_nodes, b_elems] = extract_core_mesh(boundary, 'surface');
[m_nodes, m_elems] = extract_core_mesh(membrane, 'volume');
fprintf('Assembling BEM system matrices ...'); tic;
[Lb, Mb, Lm, Mm] = laplace_2d_coll_mex(b_nodes, b_elems, m_nodes, m_elems);
fprintf('Ready in %.2f seconds.\n', toc);
%% Compute analytical and BEM solutions
m_cent = centnorm(membrane);
b_cent = centnorm(boundary);

% Position of excitation
x_exc = r0*[0.7 0.2];
ds_exc = bsxfun(@minus, b_cent(:,1:2), x_exc);
ds_exc = sqrt(dot(ds_exc, ds_exc, 2));
us_inc = -1 / (2*pi) * log(ds_exc);

dm_exc = bsxfun(@minus, m_cent(:,1:2), x_exc);
dm_exc = sqrt(dot(dm_exc, dm_exc, 2));
um_inc = -1 / (2*pi) * log(dm_exc);

us_scat = -us_inc;

ds_scat = Lb \ (Mb * us_scat);
um_scat = Mm * us_scat - Lm * ds_scat;

um_bem = um_inc + um_scat;

r_exc = sqrt(dot(x_exc, x_exc, 2));
x_exc_star = x_exc * r0^2 / r_exc^2;
dm_exc_star = bsxfun(@minus, m_cent(:,1:2), x_exc_star);
dm_exc_star = sqrt(dot(dm_exc_star, dm_exc_star, 2));

um_ana = 1/(2*pi)*(-log(dm_exc) + log(dm_exc_star) + log(r_exc/r0));
%%
err_m = norm(um_bem-um_ana) / norm(um_ana);
fprintf('Membrane error: %g\n', err_m);
%%
figure('Renderer', 'painters');
plot_mesh(membrane, um_bem);
axis equal;
caxis([0, max(um_bem)]);
figure('Renderer', 'painters');
plot_mesh(membrane, um_ana);
axis equal;
caxis([0, max(um_bem)]);