function mesh = import_bulk_mesh2(fname)
%IMPORT_BULK_MESH Import NiHu mesh from bulk file (.bdf)
%   MESH = IMPOT_BULK_MESH(FNAME)  imports the NiHu mesh from the bulk
%   file.

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

%% Parameter check
error(nargchk(1, 1, nargin, 'struct'));

%% Open and read the file
fprintf(1, 'Importing bulk file...');
file = textread(fname, '%s', 'delimiter', '\n', 'whitespace', '');
fprintf(1, 'Done\n');

%% Read node information
Nodes = zeros(0,4);
rows = strmatch('GRID*   ', file);
if ~isempty(rows)
    l = length(rows);
    fprintf(1, 'Reading %d GRID* entries...', l);
    rows = unique([rows; rows+1]);
    c = file(rows);
    % replace number1[+-]number2 to E[+1]number2
    c = regexprep(c, '(\d)([+-])(\d)', 'E$2$3');
    data = sscanf(cell2mat(c'), 'GRID*%d%16g%16g*%16g', [4, l]);
    if numel(data) ~= 4*l
        data = sscanf(cell2mat(c'), 'GRID*%d%16g%16g*%16g%*16g', [4, l]);
    end
    Nodes = [
        Nodes
        data'
        ];
    file = file(setdiff(1:length(file), rows));
    fprintf(1, 'Done\n');
end

rows = strmatch('GRID    ', file);
if ~isempty(rows)
    l = length(rows);
    fprintf(1, 'Reading %d GRID entries...', l);
    c = file(rows);
    data = sscanf(cell2mat(c'), 'GRID    %8d%8g%8g%8g', [4, l]);
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
    rows = strmatch(name, file);
    if isempty(rows)
        continue;
    end
    l = length(rows);
    fprintf(1, 'Reading %d %s entries...', l, strtrim(name));
    type = ElementTypes{net,2};
    nnod = ElementTypes{net,3};
    if nnod > 6
        rows = unique([rows; rows+1]);
    end
    ind = max(ind) + (1:l);
    c = file(rows);
    data = sscanf(cell2mat(c'), [name repmat('%8d', 1, nnod+2)], [nnod+2,l]);
    if numel(data) ~= (l*(nnod+2))
        data = sscanf(cell2mat(c'), [name repmat('%8d', 1, nnod+1)], [nnod+1,l]);
        data = [
            data(1,:)
            ones(1,l)
            data(2:end,:)
            ];
    end
    Elements(ind, [1 3 4+(1:nnod)]) = data';
    Elements(ind, 2) = type;
    Elements(ind, 4) = 1;
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