function mesh = import_bulk_mesh(fname)
%IMPORT_BULK_MESH Import NiHu mesh from bulk file (.bdf)
%   MESH = IMPOT_BULK_MESH(FNAME)  imports mesh from a Nastran Bulk file.
%
% See also: EXPORT_BULK_MESH

%   Copyright 2008-2012 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified: 2012.12.14.

%% Open and read the file
fid = fopen(fname, 'rt');
if fid == -1
    error('Cannot open bulk file %s', fname);
end
file = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
file = file{1};
fclose(fid);

%% Read GRID* information
Nodes = zeros(0,4);
name = 'GRID*   ';
rows = find(strncmp(name, file, length(name)));
if ~isempty(rows)
    l = length(rows);
    fprintf(1, 'Reading %d GRID* entries...', l);
    c1 = file(rows);
    c1 = regexprep(c1, '.{8}(.{16}).{16}(.{16})(.{1,16}).*', '$1 $2 $3 ');
    c2 = file(rows+1);
    c2 = regexprep(c2, '.{8}(.{1,16}).*', '$1');
    c = regexprep([c1 c2]', '(\d)([+-])(\d)', 'E$2$3');
    data = sscanf(cell2mat(c(:)'), '%d%g%g%g', [4, l]);
    Nodes = [
        Nodes
        data'
        ];
    file = file(setdiff(1:length(file), [rows; rows+1]));
    fprintf(1, 'Done\n');
end

%% Read GRID information
name = 'GRID    ';
rows = find(strncmp(name, file, length(name)));
if ~isempty(rows)
    l = length(rows);
    fprintf(1, 'Reading %d GRID entries...', l);
    c = file(rows);
    c = regexprep(c, '.{8}(.{8}).{8}(.{8})(.{8})(.{1,8}).*', '$1 $2 $3 $4 ');
    data = sscanf(cell2mat(c'), '%d%g%g%g', [4, l]);
    Nodes = [
        Nodes
        data'
        ];
    file = file(setdiff(1:length(file), rows));
    fprintf(1, 'Done\n');
end

%% Read element information
ElementTypes = [
    {'CROD    '}, {12}, {2}
    {'CBAR    '}, {12}, {2}
    {'CTRIA3  '}, {23}, {3}
    {'CQUAD4  '}, {24}, {4}
    {'CTETRA  '}, {34}, {4}
    {'CPENTA  '}, {36}, {6}
    {'CHEXA   '}, {38}, {8}
    ];
ind = 0;
for net = 1 : size(ElementTypes)
    name = ElementTypes{net,1};
    rows = find(strncmp(name, file, length(name)));
    if isempty(rows)
        continue;
    end
    l = length(rows);
    fprintf(1, 'Reading %d %s entries...', l, strtrim(name));
    type = ElementTypes{net,2};
    nnod = ElementTypes{net,3};
    if nnod > 6
        nnod2 = nnod - 6;
        nnod = 6;
    else
        nnod2 = 0;
    end
    ind = max(ind) + (1:l);
    c = regexprep(file(rows),...
        sprintf('.{8}(.{8}).{8}(.{%d})(.{1,8}).*', 8*(nnod-1)), '$1 $2 $3 ');
    if nnod2 ~= 0
        c2 = regexprep(file(rows+1),...
            sprintf('.{8}(.{%d})(.{1,8}).*', 8*(nnod2-1)), '$1 $2');
        c = [c c2]';
        nnod = nnod + nnod2;
        rows = [rows; rows+1]; %#ok<AGROW>
    end
    data = sscanf(cell2mat(c(:)'), repmat('%d', 1, nnod+1), [nnod+1,l]);
    Elements(ind, [1 4+(1:nnod)]) = data'; %#ok<AGROW>
    Elements(ind, 2) = type; %#ok<AGROW>
    Elements(ind, 3:4) = 1; %#ok<AGROW>
    file = file(setdiff(1:length(file), rows));
    fprintf(1, 'Done\n');
end

%% Properties and Materials
Properties = unique(Elements(:,4));
Materials = unique(Elements(:,3));

mesh.Properties = Properties;
mesh.Materials = Materials;
mesh.Nodes = Nodes;
mesh.Elements = Elements;

end