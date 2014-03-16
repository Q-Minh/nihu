clear;

%% parameters
R = 1;
H = 4;
Le = .15;
nu = .33;

%% create mesh
N = ceil(2*pi*R/Le/4);
disc = create_circle(R, N);
b = get_boundary(disc);
Nz = ceil(H/Le);
side = extrude_mesh(b, [0 0 H/Nz], Nz);
mesh = join_meshes(disc,...
    flip_elements(side),...
    translate_mesh(flip_elements(disc), [0 0 H]));
mesh = merge_coincident_nodes(mesh);
mesh = drop_unused_nodes(mesh);

%% quadratic elements and blow up
qmesh = quadratise(mesh);
[nodidx, elidx] = mesh_select(qmesh, sprintf('abs(r-%g) < 1e-2', R), 'ind');
sidenodes = qmesh.Nodes(nodidx,2:3);
r = sqrt(dot(sidenodes, sidenodes, 2));
qmesh.Nodes(nodidx,2:3) = bsxfun(@times, qmesh.Nodes(nodidx,2:3), R./r);

%% Compile mex file
mex -v CXXFLAGS="$CXXFLAGS -std=c++11 -O3 -fPIC" ...
    ../../src/tutorial/elastostatics.mex.cpp ...
    ../../src/library/lib_element.cpp ...
    ../../src/library/lib_shape.cpp ...
    ../../src/library/lib_domain.cpp ...
    -I../../src ...
    -I/usr/local/include/eigen3 ...
    -o elastostatics

%% Call NiHu to compute system matrices
[nodes, elements] = extract_core_mesh(qmesh);
tic;
[U, T] = elastostatics(nodes, elements, nu);
toc;

%% Boundary conditions
[xc, nc] = centnorm(mesh);
phi = atan2(xc(:,2), xc(:,1));
tpre = [sin(phi) -cos(phi) zeros(size(phi,1),1)];
tpre(xc(:,3) < H,:) = 0;

upind = find(xc(:,3) > 0);
updofind = union(union(3*upind-2, 3*upind-1), 3*upind);

Uc = U(updofind, updofind);
Tc = T(updofind, updofind);

ures(updofind) = (Tc - .5*eye(size(Tc))) \ (Uc * reshape(tpre(upind,:).', [], 1));
ures = reshape(ures, 3, []).';

Trans = cell2node_interp(mesh);
Trans(end,:) = Trans(end,:)/4;
ures_node = Trans * ures;

%% plot results
disp_mesh = mesh;
disp_mesh.Nodes(:,2:4) = disp_mesh.Nodes(:,2:4)-ures_node/10;
figure;
formatfig([6 10], [0 0 0 0]);
% shading interp;
quiver3(xc(:,1), xc(:,2), xc(:,3)+.5, tpre(:,1), tpre(:,2), tpre(:,3), 'LineWidth', 2, 'Color', 'black');
hold on;
plot_mesh(disp_mesh, ures_node(:,1));
view(56, 25);
light
lighting phong
axis off
print -dpng bar -r600