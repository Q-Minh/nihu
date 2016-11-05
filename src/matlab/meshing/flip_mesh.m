function mesh = flip_mesh(mesh)
%FLIP_MESH  Flip elements of a NiHu mesh
%   MESH = FLIP_ELEMENTS(MESH) flips all the elements of a NiHu mesh

%   Copyright 2008-2016 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified 2016.01.07.

mesh.Elements = flip_elements(mesh.Elements);


end % of function
