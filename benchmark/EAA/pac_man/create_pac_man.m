function [mesh, field, dirfield, dirreffield] = create_pac_man(Le)
% References:
%   [1] https://eaa-bench.mec.tuwien.ac.at/fileadmin/t/eaa/PACMAN_benchmark_case.pdf

phi0 = pi/6;
R = 1;

%% create boundary mesh
N = ceil(R * (2*pi - 2*phi0) / Le);

phi = linspace(phi0, 2*pi-phi0, N+1).';
x2 = [R * cos(phi) R * sin(phi)];
x2 = x2(1:end-1,:);

N = ceil(R / Le);

l = linspace(0, R, N+1).';
x1 = [l * cos(phi0) l * sin(phi0)];
x1 = x1(1:end-1,:);

l = linspace(R, 0, N+1).';
x3 = [l * cos(-phi0) l * sin(-phi0)];
x3 = x3(1:end-1,:);

x = [x1; x2; x3];
el = [(1:size(x,1))', [(2:size(x,1))';1]];

N = size(x,1);
mesh = create_empty_mesh();
mesh.Nodes = [(1:N)' x zeros(N,1)];

mesh.Elements = [(1:N)', repmat(ShapeSet.LinearLine.Id, N,1), ones(N,2), ...
    el];

%% create field point mesh
field = create_slab([4 4], ceil([4/Le 4/Le]));
field = translate_mesh(field, [-2 -2]);
[~, keep] = mesh_select(field, 'r > 1 | abs(phi) < pi/6', 'ind');
field.Elements = field.Elements(keep,:);
field = drop_unused_nodes(field);
field = drop_mesh_IDs(field);

%% create field point mesh for directivity
R0 = 2;
N_phi = 360*100;
dirfield = create_circle_boundary(R0, N_phi);
%dirfield = rotate_mesh(dirfield, pi/N_phi, [0 0 -1]);
xf = centnorm(dirfield);
dirfield = scale_mesh(dirfield, mean(R0./sqrt(dot(xf,xf,2))));
dirfield.Nodes(:,4) = 0;

%% Create field point mesh for reference directivity
r = 2;
Nphi = 72;
dirreffield = create_circle_boundary(r, Nphi);
%dirreffield = rotate_mesh(dirreffield, pi/N_phi, [0 0 -1]);
xf = centnorm(dirreffield);
dirreffield = scale_mesh(dirreffield, mean(r./sqrt(dot(xf,xf,2))));
dirreffield.Nodes(:,4) = 0;

end % of function create_pac_man
