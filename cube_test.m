clear;

%% geometry
elemtype = 'mixed';

model = create_brick_boundary(2, 20);
model = translate_mesh(model, [-1 -1 -1]);
if strcmp(elemtype, 'tria')
    model = quad2tria(model);
elseif strcmp(elemtype, 'mixed')
    m1 = mesh_section(model, [
        -Inf, -Inf, -Inf;
        Inf, Inf, 1e-2]);
    m2 = mesh_section(model, [
        -Inf, -Inf, -1e-2;
        Inf, Inf, Inf]);
    model = join_meshes(m1, quad2tria(m2));
    model = merge_coincident_nodes(model);
    model = drop_unused_nodes(model);
end
[c, n] = centnorm(model);
[~, ~, w, gind] = geo2gauss(model, 0);
w(gind) = w;

%% system matrices
k = 1;
[nodes, elements]  = extract_Boonen_mesh(model);
tic;
[G, H, n_eval] = Boonen13_gal_const(nodes, elements, k);
toc;
M = diag(w);

%% excitation
src = [-.2 -.3 0];
[p_inc, q_inc] = incident('point', src, c, n, k);

%% solution and error
p = (H - .5*M) \ (G*q_inc);
err = log10(abs(p-p_inc)./abs(p_inc));

%% plots
figure;
subplot(1,2,2);
plot_mesh(model, err);
shading flat;
caxis([-5 0]);
colorbar;
subplot(1,2,1);
plot_mesh(model, real(p_inc));
shading flat;
colorbar;

