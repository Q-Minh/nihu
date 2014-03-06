clear;

%% parameters
R = 1;
H = 4;
Le = .2;
nu = .33;

%% create mesh
N = ceil(2*pi*R/Le/4);
disc = create_circle(R, N);
b = get_boundary(disc);
Nz = ceil(H/Le);
side = extrude_mesh(b, [0 0 H/Nz], Nz);
mesh = join_meshes(flip_elements(disc), side, translate_mesh(disc, [0 0 H]));
mesh = merge_coincident_nodes(mesh);
mesh = drop_unused_nodes(mesh);

%% quadratic elements and blow up
qmesh = quadratise(mesh);
[nodidx, elidx] = mesh_select(qmesh, sprintf('abs(r-%g) < 1e-2', R), 'ind');
sidenodes = qmesh.Nodes(nodidx,2:3);
r = sqrt(dot(sidenodes, sidenodes, 2));
qmesh.Nodes(nodidx,2:3) = bsxfun(@times, qmesh.Nodes(nodidx,2:3), R./r);

%% Call system matrices
[nodes, elements] = extract_core_mesh(qmesh);
[L, M] = elastostatics(nodes, elements, nu);
