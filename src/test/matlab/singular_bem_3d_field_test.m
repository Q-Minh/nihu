clear;
close all;

tria = true;
faketria = ~true;

%// sphere radiator with radius R = 1, divided into 6 elements along the radius.
Lout = [1 1 1];
Lin = [.6 .6 .6];
Lop = .2;
Le = 1e-1;
mesh = create_brick(Lout, ceil(Lout/Le));
mesh = translate_mesh(mesh, -Lout/2);

eps = 1e-3;
mesh = mesh_section(mesh, [-1-eps; 1+eps] * Lin/2, 'nall');
mesh = mesh_section(mesh, [-Lop/2 -Lop/2 0; Lop/2 Lop/2 Inf]*(1+eps), 'nall');

mesh = drop_mesh_IDs(drop_unused_nodes(get_boundary(mesh)));
if (tria)
    mesh = quad2tria(mesh);
    if (faketria)
        mesh.Elements(:,2) = 124;
        mesh.Elements(:,8) = mesh.Elements(:, 7);
    end
end

field = mesh;
eps = 1e-5;
[x, n] = centnorm(mesh);
field_centers = x + eps * n;
nField = size(x, 1);
if (tria && ~faketria)
    field_elems = [32303*ones(nField, 1), repmat((0 : nField-1).', 1, 3)];
else
    field_elems = [32404*ones(nField, 1), repmat((0 : nField-1).', 1, 4)];
end

%// extract the mesh description matrices
[r_nodes, r_elem] = extract_core_mesh(mesh);

%// call C++ code at wave number k \approx pi
c = 340;
f = 38;
k = 2*pi*f/c;
[Ls, Ms, Mts, Ns, Lf, Mf] =...
    singular_bem_3d_field(r_nodes, r_elem, field_centers, field_elems, k);

Mts_num = (Lf - Ls) / eps;
diff = Mts_num - Mts;
diff = diff - diag(diag(diff));
fprintf('Diff in Mt: %f\n', ...
    norm(diff, 'fro') / norm(Mts - diag(diag(Mts)), 'fro'));

Ns_num = (Mf - Ms) / eps;
diff = Ns_num - Ns;
diff = diff - diag(diag(diff));
fprintf('Diff in N: %f\n', ...
    norm(diff, 'fro') / norm(Ns - diag(diag(Ns)), 'fro'));

%// define the incident velocity field and analytical solutions
x0 = [0 0 0];
nElem = size(mesh.Elements, 1);
qs_ana = ones(nElem, 1);

%// acoustic pressure on the surface and in the field points
ps_conv = Ms \ (Ls * qs_ana);               %// conventional
coup = 1i/k;   %// coupling constant
ps_bm = (Ms + coup * Ns) \ (Ls * qs_ana + coup * Mts * qs_ana);

ps_hs = (Ns) \ (Mts * qs_ana);

% // plot results
figure;
plot_mesh(mesh, abs(ps_conv));
figure;
plot_mesh(mesh, abs(ps_bm));
figure;
plot_mesh(mesh, abs(ps_hs));