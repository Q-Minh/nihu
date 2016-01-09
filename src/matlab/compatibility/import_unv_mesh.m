function mesh = import_unv_mesh(fname)
%IMPOR_UNV_MESH Import mesh from UNV format
%   MESH = IMPORT_UNV_MESH(FNAME) imports a mesh from unv format

% 2411: nodes
% 2412: elements

fid = fopen(fname);
if fid == -1
    error('NiHu:import_unv_mesh:invalid_file', ...
        'cannot open unv file: %s', fname);
end
data = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);
data = data{1};

blocksep = find(strcmp('-1', data));
blockind = [blocksep(1:2:end)+1 blocksep(2:2:end)-1];
blockid = str2double(data(blockind(:,1)));

nodeblock = find(blockid == 2411);
elemblock = find(blockid == 2412);
unusedblocks = setdiff(blockid, [2411 2412]);
if ~isempty(unusedblocks)
    warning('NiHu:read_unv_mesh', ...
        sprintf('Block not imported: %g\n', unusedblocks)); %#ok<SPWRN>
end

nodedef = data(blockind(nodeblock,1)+1 : blockind(nodeblock,2));
nodedef = strrep(nodedef, 'D', 'E');
elemdef = data(blockind(elemblock,1)+1 : blockind(elemblock,2));

% read nodes
nNodes = size(nodedef,1)/2;
nodid = zeros(nNodes,1);
coord = zeros(nNodes,3);
i = 0;
line = 1;
while line < length(nodedef)
    [res1, count1] = sscanf(nodedef{line}, '%u', [4 1]);
    if count1 ~= 4
        break;
    end
    [res2, count2] = sscanf(nodedef{line+1}, '%g', [3 1]);
    if count2 ~= 3
        break;
    end
    i = i + 1;
    nodid(i) = res1(1);
    coord(i,:) = res2;
    line = line + 2;
end

% read elements

elems = zeros(0,3);
elemids = zeros(0,1);
elemtypes = zeros(0,1);
i = 0;
line = 1;
while line < length(elemdef)
    [res1, count1] = sscanf(elemdef{line}, '%u', [6 1]);
    if count1 ~= 6
        break;
    end
    nnodes = res1(6);
    [res2, count2] = sscanf(elemdef{line+1}, '%u', [nnodes 1]);
    if count2 ~= nnodes
        break;
    end
    i = i + 1;
    elemids(i,1) = res1(1);
    elemtypes(i,1) = res1(2);
    elems(i,1:nnodes) = res2;
    line = line + 2;
end

mesh = create_empty_mesh();
mesh.Nodes = [nodid(:) coord];
mesh.Elements = zeros(size(elems,1), 7);
mesh.Elements(:,4+(1:size(elems,2))) = elems;
tri = sum(mesh.Elements ~= 0, 2) == 3;
qua = sum(mesh.Elements ~= 0, 2) == 4;
qtri = elemtypes == 92;
mesh.Elements(tri,2) = ShapeSet.LinearTria.Id;
mesh.Elements(qua,2) = ShapeSet.LinearQuad.Id;
mesh.Elements(qtri,2) = ShapeSet.QuadraticTria.Id;
mesh.Elements(:,3:4) = 1;
mesh.Elements(:,1) = elemids;

end % of function