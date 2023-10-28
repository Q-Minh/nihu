% circular_membrane.m
% Circular membrane under static load at a single point

% Comparison of boundary element method and Green's function obtained using
% the method of mirror images for the circular membrane

clear;

r0 = 0.4;
n0 = 6;
use_quadratic_elems = true;

%% Create the membrane and boundary meshes
membrane = create_circle(r0, n0);
boundary = drop_mesh_IDs(drop_unused_nodes(get_boundary(membrane)));
% Interior problem, normals are flipped
boundary = flip_mesh(boundary);
if (use_quadratic_elems)
    boundary = quadratise(boundary);
    xyz = boundary.Nodes(:,2:4);
    boundary.Nodes(:,2:4) = bsxfun(@times, xyz, r0./sqrt(dot(xyz, xyz, 2)));
end

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
x_exc = r0*[0.3 0.2];
db_exc = bsxfun(@minus, b_cent(:,1:2), x_exc);
db_exc = sqrt(dot(db_exc, db_exc, 2));
ub_inc = -1 / (2*pi) * log(db_exc);

dm_exc = bsxfun(@minus, m_cent(:,1:2), x_exc);
dm_exc = sqrt(dot(dm_exc, dm_exc, 2));
um_inc = -1 / (2*pi) * log(dm_exc);

% Scattered field on the boundary
ub_scat = -ub_inc;

fprintf('Solving BEM system ... '); tic;
qb_scat = Lb \ (Mb * ub_scat);
fprintf('Ready in %.2f seconds.\n', toc);

fprintf('Computing BEM field ... '); tic;
um_scat = Mm * ub_scat - Lm * qb_scat;
fprintf('Ready in %.2f seconds.\n', toc);

% Total field on the membrane computed by the BEM
um_bem = um_inc + um_scat;

%% Analytical solution using mirror images
r_exc = sqrt(dot(x_exc, x_exc, 2));
x_exc_star = x_exc * r0^2 / r_exc^2;
dm_exc_star = bsxfun(@minus, m_cent(:,1:2), x_exc_star);
dm_exc_star = sqrt(dot(dm_exc_star, dm_exc_star, 2));

um_ana = 1/(2*pi)*(-log(dm_exc) + log(dm_exc_star) + log(r_exc/r0));
%% Compute error compared to analytical solution
err_m = norm(um_bem-um_ana) / norm(um_ana);
fprintf('Relative error of displacement on membrane: %g\n', err_m);
%% Create the figure
fig = figure;
scale = .3;
plot_mesh(membrane, um_bem);
hold on;
h_exc = plot3(x_exc(1), x_exc(2), 0, 'r+', 'LineWidth', 2);
axis equal;
caxis([0, max(um_bem)]);
axis off;
cb = colorbar('Location', 'EastOutside');
ylabel(cb,'Displacement [-]', 'FontSize', 10);
title('Total displacement on the membrane');
%print(fig, '-dpng', '-r100', 'bench_circular_membrane.png');