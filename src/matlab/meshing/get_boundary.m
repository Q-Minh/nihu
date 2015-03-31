function [boundary, elemind] = get_boundary(model)
%GET_BOUNDARY  Extract domain boundary
%   [BOUNDARY, PARENTS] = GET_BOUNDARY extracts and returns the boundary
%   elements of the NiHu mesh MESH. The boundary of the TRIA and QUAD
%   elements are LINEs. The boundary of TETRA elements are TRIA, etc.
%   The mesh BOUNDARY contains all the nodes of the initial mesh.
%   The elements of BOUNDARY have the same property and material as their
%   parent elements. The output argument vector PARENTS contains the
%   indices of each boundary element's parent in the matrix MESH.Elements.
%
% Example:
%   mesh = create_brick(1, 10);
%   x = mesh.Nodes(:,2);
%   x = mean(x(mesh.Elements(:,5:8)), 2);
%   [bou, parent] = get_boundary(mesh);
%   plot_mesh(bou, x(parent));
%
% See also: GET_FACES, GET_FREE_FACES

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified 2012.12.19.

% Extract boundary
boundary.Materials = model.Materials;
boundary.Properties = model.Properties;
boundary.Nodes = model.Nodes;
faces = get_faces(model.Elements);
faces = get_free_faces(faces);
boundary.Elements = zeros(size(faces,1),0);
boundary.Elements(:,[2 4+(1:size(faces,2)-2)]) = faces(:,2:end);
[~, elemind] = ismember(faces(:,1), model.Elements(:,1));
if ~isempty(boundary.Elements)
    boundary.Elements(:,[3 4]) = model.Elements(elemind,[3 4]);
    boundary.Elements(:,1) = 1:size(boundary.Elements,1);
end

end
