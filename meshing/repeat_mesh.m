function mesh2 = repeat_mesh(mesh, dir, nrep)
%REPEAT_MESH Repeat (Copy) mesh along a given vector
%   MESH = REPEAT(MESH, DIR, NREP) repeats the mesh MESH
%   NREP times so that each new repetition is the original mesh translated
%   along the direction vector DIR.
%
% See also: TRANSLATE_MESH, SCALE_MESH, ROTATE_MESH, EXTRUDE_MESH,
% REVOLVE_MESH, REFLECT_MESH

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 02.12.2009


%
mesh2 = create_empty_mesh();
for iRep = 0 : nrep
    mesh2 = join_meshes(mesh2, translate_mesh(mesh, iRep * dir));
end

end

