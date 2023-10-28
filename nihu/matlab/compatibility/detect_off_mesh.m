function is_off = detect_off_mesh(filename)
%DETECT_OFF_MESH Detect if mesh is in OFF format
%   IS_OFF = DETECT_OFF_MESH(FILENAME)
%
% See also: import_off_mesh export_off_mesh

%   Copyright 2008-2013 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Networked Systems and Services

% Last modified: 2013.11.08.

fid = fopen(filename, 'rt');
if (fid == -1)
    error('NiHu:file_open', 'Cannot open file %s', filename);
end
s = fgetl(fid);
is_off = strcmp(s, 'OFF');

fclose(fid);

end
