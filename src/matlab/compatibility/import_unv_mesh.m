function mesh = import_unv_mesh(fname)
%IMPORT_UNV_MESH Import mesh from UNV format
%   MESH = IMPORT_UNV_MESH(FNAME) imports a mesh from unv format

% 163: units
% 2411: nodes
% 2412: elements
fid = fopen(fname);
if fid == -1
    error('nihu:runtime_error', 'Could not open unv file ''%s'' for reading', fname);
end
data = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);
data = data{1};

blocksep = find(strcmp('-1', data));
blockind = [blocksep(1:2:end)+1 blocksep(2:2:end)-1];
blockid = str2double(data(blockind(:,1)));

unitblock = find(blockid == 164);
nodeblock = find(blockid == 2411);
elemblock = find(blockid == 2412);

nodedef = data(blockind(nodeblock,1)+1 : blockind(nodeblock,2));
nodedef = strrep(nodedef, 'D', 'E');
elemdef = data(blockind(elemblock,1)+1 : blockind(elemblock,2));

% read units
if ~isempty(unitblock)
    unitblock = unitblock(end);
    units = process_units(data(blockind(unitblock,1)+1 : blockind(unitblock,2)));
else
    units = struct('length_ratio', 1, 'force_ratio', 1, 'temp_ratio', 1);
end

% read nodes
nNode = length(nodedef)/2;
nodid = zeros(nNode,1);
coord = zeros(nNode,3);
for i = 1 : nNode
    line = 2*i-1;
    [res1, count1] = sscanf(nodedef{line}, '%u', [4 1]);
    if count1 ~= 4
        break;
    end
    [res2, count2] = sscanf(nodedef{line+1}, '%g', [3 1]);
    if count2 ~= 3
        break;
    end
    nodid(i) = res1(1);
    coord(i,:) = res2;
end

% apply units to coordinates
coord = coord / units.length_ratio;

% read elements
nElem = length(elemdef)/2;
elems = zeros(nElem,3);
elemids = zeros(nElem,1);
elemtypes = zeros(nElem,1);
for i = 1 : nElem
    line = 2*i-1;
    [res1, count1] = sscanf(elemdef{line}, '%u', [6 1]);
    if count1 ~= 6
        break;
    end
    nnodes = res1(6);
    [res2, count2] = sscanf(elemdef{line+1}, '%u', [nnodes 1]);
    if count2 ~= nnodes
        break;
    end
    elemids(i,1) = res1(1);
    elemtypes(i,1) = res1(2);
    elems(i,1:nnodes) = res2;
end

mesh = create_empty_mesh();
mesh.Nodes = [nodid(:) coord];
mesh.Elements = zeros(size(elems,1), 7);
mesh.Elements(:,4+(1:size(elems,2))) = elems;
lin = sum(mesh.Elements ~= 0, 2) == 2;
tri = sum(mesh.Elements ~= 0, 2) == 3;
qua = sum(mesh.Elements ~= 0, 2) == 4;
qtri = elemtypes == 92;
mesh.Elements(lin,2) = ShapeSet.LinearLine.Id;
mesh.Elements(tri,2) = ShapeSet.LinearTria.Id;
mesh.Elements(qua,2) = ShapeSet.LinearQuad.Id;
mesh.Elements(qtri,2) = ShapeSet.QuadraticTria.Id;
mesh.Elements(:,3:4) = 1;
mesh.Elements(:,1) = elemids;

end % of function


function units = process_units(data)
code = sscanf(data{1}, '%u');
ratios = sscanf(strrep(data{2}, 'D', 'E'), '%f');
units.length_ratio = ratios(1);
units.force_ratio = ratios(2);
units.temp_ratio = ratios(3);
units.temp_offset = sscanf(strrep(data{3}, 'D', 'E'), '%f');
end
