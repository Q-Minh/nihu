function is_bulk = detect_bulk_mesh(fname)
%DETECT_BULK_MESH Detect if mesh is in bulk format
%   IS_BULK = DETECT_BULK_MESH(FNAME)  detects if FILENAME is a BULK mesh
%
% See also: EXPORT_BULK_MESH IMPORT_BULK_MESH

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

pattern = 'GRID';
is_bulk = any(strncmp(pattern, ftext, length(pattern)));

end
