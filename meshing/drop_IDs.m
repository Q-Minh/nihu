function Elements = drop_IDs(model)
%DROP_IDS   Get rid of material, property, element and node IDs
%   ELEMENTS = DROP_IDS(MESH) returns an updated version of the matrix
%   MODEL.Elements. This new matrix refers NOT to the material, property
%   and node IDs, but the corresponding indices in the matrices
%   MODEL.Nodeds, MODEL.Properties and MODEL.Materials.

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last updated: 02.12.2009.

%% Parameter check
if size(model.Elements,1) == 0
    Elements = model.Elements;
    return
end

%% Create permutation Index vectors
% for nodes
NodeID = model.Nodes(:,1);
NodeIndex(NodeID+1) = (1:size(model.Nodes,1))+1;
NodeIndex(1) = 1;
% for properties
PropertyID = model.Properties(:,1);
PropertyIndex(PropertyID+1) = (1:size(model.Properties,1))+1;
PropertyIndex(1) = 1;
% for materials
MaterialID = model.Materials(:,1);
MaterialIndex(MaterialID+1) = (1:size(model.Materials,1))+1;
MaterialIndex(1) = 1;

%% Update matrix model.Elements using the permutation vectors
Elements = model.Elements;
Elements(:,1) = 1:size(Elements,1);
Elements(:,3) = MaterialIndex(Elements(:,3)+1)-1;
Elements(:,4) = PropertyIndex(Elements(:,4)+1)-1;
Elements(:,5:end) = NodeIndex(Elements(:,5:end)+1)-1;
end
