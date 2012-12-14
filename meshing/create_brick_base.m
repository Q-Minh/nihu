function brick = create_brick_base(N)
%CREATE_BRICK_BASE Create a basic brick mesh. (NiHu / meshing)
% Creates a basic brick mesh which is used for further
% transformations. The basic brick is hexahedron with edge length 2 
% centered at (0,0,0) and having N(1), N(2) and N(3) elements along the
% edges.
%
% See also: CREATE_BRICK, CREATE_SLAB_BASE

brick = translate_mesh(...
    extrude_mesh(create_slab_base(N(1:2)), [0 0 2/N(3)], N(3)), [0 0 -1]);

end % create_brick_base
