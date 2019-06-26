function [mesh, field, k, q_surf] = create_pac_man_validation(Le, freq)

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
r = 2;
Nphi = 72;
field = rotate_mesh(create_circle_boundary(r, Nphi), pi/Nphi, [0 0 -1]);
xf = centnorm(field);
field = scale_mesh(field, mean(2./sqrt(dot(xf,xf,2))));
field.Nodes(:,4) = 0;

%% create surface velocity excitation
rho = 1.2041;
c = 343.21;
V0 = 0.1;

% freq = 4000;
om = 2* pi * freq;
k = om / c;

v_surf = zeros(N,1);
[~, sel] = mesh_select(mesh, 'abs(r-1) < 1e-4', 'ind');
v_surf(sel) = V0;

q_surf = -1i*om*rho*v_surf;

end
