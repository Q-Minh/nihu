function [mesh, points] = create_helmholtz_bem(Le, trias)

if nargin < 2
    trias = false;
end

L = [1.0 1.0 1.0];
N = ceil(L/Le);
mesh = create_brick(L, N);
mesh = translate_mesh(mesh, -L/2);

mesh0 = mesh;

bricks = {
    [-0.3, 0.3], [-0.3, 0.3], [-0.3, 0.3];
    [-0.1, 0.1], [-0.1, 0.1], [0.3, 0.5]
    };

for b = 1 : size(bricks,1)
    limits = [bricks{b,1}' bricks{b,2}' bricks{b,3}'];
    limits = bsxfun(@plus, limits, Le/10*[-1; 1]);
    selmode = 'nall';
    mesh = mesh_section(mesh, limits, selmode);
end

mesh = drop_mesh_IDs(mesh);

mesh = drop_mesh_IDs(drop_unused_nodes(get_boundary(mesh)));
if (trias)
    mesh = quad2tria(mesh);
end

points = [0, 0, L(3)/2 + .2];

end % of function create_helmholtz_bem
