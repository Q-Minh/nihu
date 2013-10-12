function mesh = import_off_mesh(filename)
%IMPORT_OFF_MESH import NiHu mesh from OFF format
%   MESH = IMPORT_OFF_MESH(FILENAME) imports the NiHu mesh from OFF format.
%
% See also: import_gmsh_mesh import_bulk_mesh

%   Copyright 2008-2013 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Networked Systems and Services

% Last modified: 2013.10.11.

fid = fopen(filename, 'rt');
if (fid == -1)
    error('NiHu:file_open', 'Cannot open file %s', filename);
end
s = fgetl(fid);
if ~strcmp(s, 'OFF')
    error('NiHu:invalid_data', 'invalid OFF file');
end
data = fscanf(fid, '%d', 3);
nnodes = data(1);
nelements = data(2);

% read nodes
nodes = fscanf(fid, '%g', nnodes*3);
nodes = reshape(nodes, 3, []).';

% read elements
elements = zeros(nelements,0);
for e = 1 : nelements
	n = fscanf(fid, '%u', 1);
	el = fscanf(fid, '%u', n);
	elements(e,1:(1+n)) = [n el(:)'+1];
end
fclose(fid);

% convert to NiHu format
mesh.Nodes(:,2:4) = nodes;
mesh.Nodes(:,1) = 1:nnodes;
mesh.Elements(:,5:7) = elements(:,2:4);
mesh.Elements(:,1) = 1:nelements;
mesh.Elements(:,2) = 23;
mesh.Elements(:,3:4) = 1;
[mesh.Materials, mesh.Properties] = default_mat_prop();

end

