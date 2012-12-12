function slab = create_slab_base(N)
% CREATE_SLAB_BASE Creates basic slab mesh. (NiHu / meshing)
    slab = translate_mesh(...
        extrude_mesh(create_line_base(N(1)), [0 2/N(2) 0], N(2)), [0 -1 0]);
end % create_slab_base