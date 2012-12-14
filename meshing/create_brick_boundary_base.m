function mesh = create_brick_boundary_base(N)
% CREATE_BRICK_BOUNDARY_BASE Creates a basic brick mesh. (NiHu / meshing)
% Creates a basic brick mesh which is used for further
% transformations. The basic brick is hexahedron with edge length 2
% centered at (0,0,0) and having N(1), N(2) and N(3) elements along the
% edges.
%
% See also: CREATE_BRICK_BOUNDARY, CREATE_SLAB_BASE

corners = [
    -1 -1 -1
    +1 -1 -1
    +1 +1 -1
    -1 +1 -1
    -1 -1 +1
    +1 -1 +1
    +1 +1 +1
    -1 +1 +1
    ];

faces = [
    1 4 3 2
    5 6 7 8
    1 2 6 5
    2 3 7 6
    3 4 8 7
    4 1 5 8
    ];

divs = [
    2 1
    1 2
    1 3
    2 3
    1 3
    2 3
    ];

mesh = create_empty_mesh();
for i = 1 : size(faces,1)
    mesh = join_meshes(mesh,...
        create_slab(corners(faces(i,:),:), N(divs(i,:))));
end

mesh = drop_mesh_IDs(merge_coincident_nodes(mesh, min(1./(10*N))));

end % create_brick_boundary_base
