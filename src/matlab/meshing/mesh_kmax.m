function k = mesh_kmax(model, ratio)
%MESH_KMAX Maximal wave numbers of a NiHu mesh (NiHu / meshing)
%   k = mesh_kmax(model, ratio) returns the maximal wave number for each
%   element of the mesh given by model.
% Parameters:
%   model : mesh structure
%   ratio : number of elements per wavelength (usually set to 6-8)
%   k     : nEx1 vector with wave numbers for each element

%   Copyright 2009-2012 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last updated: 2012.12.19.

% default ratio
if nargin < 2
    ratio = 8;
end

d = mesh_edge_size(model);
k = (2*pi) ./ (d * ratio);

end
