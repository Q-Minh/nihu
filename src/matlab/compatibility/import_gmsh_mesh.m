function [mesh, info] = import_gmsh_mesh(fname)
%IMPORT_GMSH_MESH Import NiHu mesh from gmsh file
%   [MESH, INFO] = IMPOT_GMSH_MESH(FNAME)  imports mesh from a gmsh file
%
% See also: EXPORT_GMSH_MESH IMPORT_BULK_MESH EXPORT_BULK_MESH

%   Copyright 2008-2013 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Networked Systems and Services

% Last modified: 2013.10.11.

fid = fopen(fname);
if fid == -1
    error('NiHu:file_open', 'Cannot open bulk file %s', fname);
end
ftext = textscan(fid, '%s', 'Delimiter', '\n');
ftext = ftext{1};
fclose(fid);

% Commands start with $ sign
rxp_command = '^\$';

r = regexp(ftext, rxp_command, 'start');

k = 1;
commands = cell(0,2);


for i = 1 : numel(r)
    if ~isempty(r{i})
        commands{k, 1} = ftext{i};
        commands{k, 2} = i;
        k = k+1;
    end
end

% Obtain info
info = read_info(ftext, commands);

% Read nodes
mesh.Nodes = read_nodes(ftext, commands);

% Read elements
mesh.Elements = read_elements(ftext, commands);

% Create default materials and properties
[mesh.Materials, mesh.Properties] = default_mat_prop();


end

function info = read_info(ftext, commands)
    start_exp = '^\$MeshFormat';
    %stop_exp = '^\$EndMeshFormat';
    
    start_ind = regexp(commands(:,1), start_exp, 'start', 'once');
    start_ind = find(~cellfun(@isempty, start_ind), 1, 'first');
    
    % If no info found
    if isempty(start_ind)
        info = [];
        return
    end
    
    % Starting line of information
    start_line = ftext{commands{start_ind,2}+1};
    
    % Process the structure
    data = sscanf(start_line, '%f', 3);
    
    info.version = data(1);
    if data(2) == 0
        info.file_type = 'ASCII';
    else
        info.file_type = data(2);
    end
    info.data_size = data(3);
end

function nodes = read_nodes(ftext, commands)
    start_exp = '^\$Nodes';
    stop_exp = '^\$EndNodes';
    
    start_ind = regexp(commands(:,1), start_exp, 'start', 'once');
    start_ind = find(~cellfun(@isempty, start_ind), 1, 'first');
    
    stop_ind = regexp(commands(:,1), stop_exp, 'start', 'once');
    stop_ind = find(~cellfun(@isempty, stop_ind), 1, 'first');
    
    % If no info found
    if isempty(start_ind) || isempty(stop_ind)
        nodes = [];
        return
    end
    
    start_line_ind = commands{start_ind,2}+2;
    stop_line_ind = commands{stop_ind,2}-1;
    
    nNodes = stop_line_ind - start_line_ind + 1;
    nodes = zeros(nNodes, 4);
    
    % ElemID, x, y, z format
    for iNode = 1 : nNodes
        nodes(iNode, :) = sscanf(ftext{start_line_ind + iNode -1}, '%f', 4);
    end
end

function elements = read_elements(ftext, commands)
    start_exp = '^\$Elements';
    stop_exp = '^\$EndElements';
    
    start_ind = regexp(commands(:,1), start_exp, 'start', 'once');
    start_ind = find(~cellfun(@isempty, start_ind), 1, 'first');
    
    stop_ind = regexp(commands(:,1), stop_exp, 'start', 'once');
    stop_ind = find(~cellfun(@isempty, stop_ind), 1, 'first');
    
    % If no info found
    if isempty(start_ind) || isempty(stop_ind)
        elements = [];
        return
    end
    
    start_line_ind = commands{start_ind,2}+2;
    stop_line_ind = commands{stop_ind,2}-1;
    
    nElem = stop_line_ind - start_line_ind + 1;
    elements = zeros(nElem, 12);
    jElem = 0;
    for iElem = 1 : nElem
        data = sscanf(ftext{start_line_ind + iElem - 1}, '%f', Inf);
        if (numel(data) < 2)
            continue
        end
        % Select elem type
        switch data(2)
            % 3-node tria
            case 2
                elements(jElem+1, 1:4) = [data(1) 23 1 1];
                elements(jElem+1, 4+(1:3)) = data(end-2 : end);
                jElem = jElem+1;
            % 4-node quadrangle
            % TODO: check node ordering 
            case 3
                elements(jElem+1, 1:4) = [data(1) 24 1 1];
                elements(jElem+1, 4+(1:4)) = data(end-3 : end);
                jElem = jElem+1;
            % 4-node tetrahedron
            case 4
                elements(jElem+1, 1:4) = [data(1) 34 1 1];
                elements(jElem+1, 4+(1:4)) = data(end-3 : end);
            % 8-node hexahedron 
            % TODO: check node ordering
            case 5
                elements(jElem+1, 1:4) = [data(1) 38 1 1];
                elements(jElem+1, 4+(1:8)) = data(end-7 : end);
            % 6-node prism
            % TODO: check node ordering
            case 6
                elements(jElem+1, 1:4) = [data(1) 36 1 1];
                elements(jElem+1, 4+(1:6)) = data(end-5 : end);
        end
    end
    % Drop unread elements
    elements = elements(1 : jElem, :);
end

