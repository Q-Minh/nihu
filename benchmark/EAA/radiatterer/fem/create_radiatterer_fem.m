function [mesh, points] = create_radiatterer_fem(Le, Ladd)

L = [2.5 2.0 1.7];
N = ceil((L + 2*Ladd)/Le);
mesh = create_brick(L + 2*Ladd, N);
mesh = translate_mesh(mesh, -Ladd*[1 1 1]);

mesh0 = mesh;

bricks = {
    [0.2, 1.4], [0.0, 1.2], [0.4, 1.2]
    [2.3, 2.5], [0.2, 0.4], [0.5, 0.7]
    [2.3, 2.5], [0.7, 0.8], [1.0, 1.1]
    [1.6, 2.3], [0.1, 1.2], [0.2, 1.5]
    [0.2, 1.0], [1.4, 1.8], [0.0, 1.7]
    [1.2, 1.4], [1.4, 2.0], [1.6, 1.7]
    [1.6, 1.8], [1.4, 2.0], [1.6, 1.7]
    [2.0, 2.2], [1.4, 2.0], [1.6, 1.7]
    [2.4, 2.5], [1.4, 2.0], [1.6, 1.7]
    [1.2, 1.4], [1.4, 2.0], [1.2, 1.4]
    [1.6, 1.8], [1.4, 2.0], [1.2, 1.4]
    [2.0, 2.2], [1.4, 2.0], [1.2, 1.4]
    [2.4, 2.5], [1.4, 2.0], [1.2, 1.4]
    [1.2, 1.4], [1.4, 2.0], [0.6, 1.0]
    [1.6, 1.8], [1.4, 2.0], [0.6, 1.0]
    [2.0, 2.2], [1.4, 2.0], [0.6, 1.0]
    [2.4, 2.5], [1.4, 2.0], [0.6, 1.0]
    [1.2, 1.4], [1.4, 2.0], [0.0, 0.4]
    [1.6, 1.8], [1.4, 2.0], [0.0, 0.4]
    [2.0, 2.2], [1.4, 2.0], [0.0, 0.4]
    [2.4, 2.5], [1.4, 2.0], [0.0, 0.4]
    };

limits = [0 0 0; 2.5 2.0 1.7];
limits = bsxfun(@plus, limits, Le/10*[-1; 1]);
mesh = drop_unused_nodes(mesh_section(mesh, limits, 'nall'));

s0 = size(mesh.Elements,1);
for b = 1 : length(bricks)
    
    
    limits = [bricks{b,1}' bricks{b,2}' bricks{b,3}'];
    limits = bsxfun(@plus, limits, Le/10*[-1; 1]);
    selmode = 'all';
    mesh_add = drop_unused_nodes(mesh_section(mesh0, limits, selmode));
    mesh = merge_coincident_nodes(join_meshes(mesh, mesh_add), Le/10);
    
end

mesh = drop_mesh_IDs(mesh);

points = [
    -0.3  0.3  0.0
     0.2 -0.3  0.4
     0.0  0.0 -0.3
     0.5  1.6  0.8
     0.6  0.5  0.8
     2.0  0.5  0.8
    ];

end
