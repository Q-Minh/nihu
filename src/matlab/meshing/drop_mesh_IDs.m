function mesh = drop_mesh_IDs(mesh)
%DROP_MESH_IDS   Get rid of material, property, element and node IDs
%   MESH = DROP_MESH_IDS(MESH) returns an updated version of the mesh
%   MESH. This new matrix refers NOT to the material, property
%   and node IDs, but the corresponding indices in the matrices
%   MESH.Nodeds, MESH.Properties and MESH.Materials.

%   Copyright 2008-2012 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last updated: 2012.12.14.

mesh.Elements = drop_IDs(mesh);
mesh.Nodes(:,1) = 1 : size(mesh.Nodes,1);
mesh.Properties(:,1) = 1 : size(mesh.Properties,1);
mesh.Materials(:,1) = 1 : size(mesh.Materials,1);

end

