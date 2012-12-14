function mesh = create_slab_boundary_base(N)
% CREATE_SLAB_BOUNDARY_BASE Creates a basic slab mesh. (NiHu / meshing)
% Creates a basic slab mesh which is used for further
% transformations. The basic slab is hexahedron with edge length 2
% centered at (0,0,0) and having N(1), N(2) and N(3) elements along the
% edges.
%
% See also: CREATE_BRICK_BOUNDARY, CREATE_SLAB_BASE

corners = [
    -1 -1
    +1 -1
    +1 +1
    -1 +1
    ];

faces = [
    1 2
    2 3
    3 4
    4 1
    ];

divs = [1 2 1 2 ];

mesh = create_empty_mesh();
for i = 1 : size(faces,1)
    mesh = join_meshes(mesh,...
        create_line(corners(faces(i,:),:), N(divs(i))));
end

mesh = drop_mesh_IDs(merge_coincident_nodes(mesh, min(1./(10*N))));

end % create_slab_boundary_base
