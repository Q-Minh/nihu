function mesh = create_empty_mesh
%CREATE_EMPTY_MESH Create an empty NiHu mesh (NiHu / meshing)

mesh.Nodes = zeros(0,4);
mesh.Elements = zeros(0,6);
[mesh.Materials, mesh.Properties] = default_mat_prop();

end

