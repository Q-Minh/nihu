function [mesh, points] = create_radiatterer(Le)

L = [2.5 2.0 1.7];
N = ceil(L/Le);
mesh = create_brick(L, N);

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

s0 = size(mesh.Elements,1);
for b = 1 : length(bricks)
    limits = [bricks{b,1}' bricks{b,2}' bricks{b,3}'];
    limits = bsxfun(@plus, limits, Le/10*[-1; 1]);
    selmode = 'nall';
    mesh = mesh_section(mesh, limits, selmode);
    s = size(mesh.Elements,1);
    rem = s0 - s;
    rem2 = diff(bricks{b,1}) * diff(bricks{b,2}) * diff(bricks{b,3}) / Le^3;
    if (round(rem) ~= round(rem2))
        error('NiHu');
    end
    s0 = s;
end

mesh = get_boundary(mesh);
mesh = drop_unused_nodes(merge_coincident_nodes(mesh));
mesh = drop_mesh_IDs(mesh);

points = [
    -.3 .3 0
    .2 -.3 .3
    0 0 -.3
    .5 1.6 .8
    .6 .5 .8
    2 .5 .8
    ];

end
