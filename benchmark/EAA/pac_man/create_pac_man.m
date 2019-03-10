function [mesh, field, xct] = create_pac_man(Le)
phi0 = pi/6;
R = 1;
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

xct = zeros(N,1);
[~, sel] = mesh_select(mesh, 'abs(r-1) < 1e-4', 'ind');
xct(sel) = 1;

field = create_slab([4 4], ceil([4/Le 4/Le]));
field = translate_mesh(field, [-2 -2]);
[~, keep] = mesh_select(field, 'r > 1 | abs(phi) < pi/6', 'ind');
field.Elements = field.Elements(keep,:);
field = drop_unused_nodes(field);
field = drop_mesh_IDs(field);

end
