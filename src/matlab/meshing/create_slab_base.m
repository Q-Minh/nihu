function slab = create_slab_base(N)
%CREATE_SLAB_BASE Create basic slab mesh. (NiHu / meshing)

l = create_line_base(N(1));
slab = translate_mesh(extrude_mesh(l, [0 2/N(2) 0], N(2)), [0 -1 0]);
end % create_slab_base
