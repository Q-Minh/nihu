function is_gmsh = detect_gmsh_mesh(fname)
%DETECT_GMSH_MESH Detect if mesh is in gmsh format
%   IS_GMSH = DETECT_GMSH_MESH(FNAME)  detects if FILENAME is a GMSH mesh
%
% See also: EXPORT_GMSH_MESH IMPORT_GMSH_MESH

%   Copyright 2008-2013 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Networked Systems and Services

% Last modified: 2013.11.08.

fid = fopen(fname);
if fid == -1
    error('NiHu:file_open', 'Cannot open bulk file %s', fname);
end
ftext = textscan(fid, '%s', 'Delimiter', '\n');
ftext = ftext{1};
fclose(fid);

% Commands start with $ sign
is_gmsh = any(strcmp('$MeshFormat', ftext));

end

