function mesh = create_empty_mesh
%CREATE_EMPTY_MESH Create an empty NiHu mesh (NiHu / meshing)

mesh.Nodes = zeros(0,4);
mesh.Elements = zeros(0,6);
mesh.Materials = [1 1 1 1 0 0];
mesh.Properties = [1 1 0 0 0 0];

end

