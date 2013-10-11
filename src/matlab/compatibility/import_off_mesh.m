function mesh = import_off_mesh(filename)
%IMPORT_OFF_MESH import NiHu mesh from OFF format
%   MESH = IMPORT_OFF_MESH(FILENAME) imports the NiHu mesh from OFF format.
%
% See also: import_gmsh_mesh import_bulk_mesh
%
% TODO: does not handle quad elements, only tria

fid = fopen(filename, 'rt');
if (fid == -1)
    error('canot open file %s', filename);
end
s = fgetl(fid);
if ~strcmp(s, 'OFF')
    error('invalid OFF file');
end
data = fscanf(fid, '%d', 3);
nnodes = data(1);
nelements = data(2);

nodes = fscanf(fid, '%g', nnodes*3);
nodes = reshape(nodes, 3, [])';
elements = fscanf(fid, '%u', nelements*4);
elements = reshape(elements, 4, [])';
elements(:,2:end) = elements(:,2:end)+1;
fclose(fid);

mesh.Nodes(:,2:4) = nodes;
mesh.Nodes(:,1) = 1:nnodes;
mesh.Elements(:,5:7) = elements(:,2:4);
mesh.Elements(:,1) = 1:nelements;
mesh.Elements(:,2) = 23;
mesh.Elements(:,3:4) = 1;
[mesh.Materials, mesh.Properties] = default_mat_prop();

end