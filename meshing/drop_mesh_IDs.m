function mesh = drop_mesh_IDs(mesh)

mesh.Elements = drop_IDs(mesh);
mesh.Nodes(:,1) = 1 : size(mesh.Nodes,1);
mesh.Properties(:,1) = 1 : size(mesh.Properties,1);
mesh.Materials(:,1) = 1 : size(mesh.Materials,1);

end

