function export_gmsh_mesh(mesh, info, fname)
%EXPORT_GMSH_MESH Export NiHu mesh to gmsh format
%
% See also: EXPORT_BULK_MESH IMPORT_BULK_MESH IMPORT_GMSH_MESH

%   Copyright 2008-2013 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified: 2013.09.22.


elements = drop_IDs(mesh);
nodes = [(1 : size(mesh.Nodes, 1)).', mesh.Nodes(:, 2:4)];

% Create default info if needed
if (isempty(info))
    info.version = 2.2;             % Version as of 2013-07-02
    info.file_type = 'ASCII';       % Use ASCII text
    info.data_size = 8;             % Use 8-byte doubles
end

% Create file and fill with content
fid = fopen(fname, 'w');            % Open file
write_info(fid, info);              % Write mesh format info
write_nodes(fid, nodes);            % Write the nodes
write_elements(fid, elements);      % Write the nodes
fclose(fid);                        % Close file

end % of export_gmsh_mesh function


function write_info(fid, info)
    fprintf(fid, '$MeshFormat\n');
    % Convert type name into value
    switch info.file_type
        case 'ASCII'
            type = 0;
        otherwise
            type = info.file_type;
    end
    fprintf(fid, '%g %d %d\n', info.version, type, info.data_size);
    fprintf(fid, '$EndMeshFormat\n');
end

function write_nodes(fid, nodes)
    fprintf(fid, '$Nodes\n');
    
    nNodes = size(nodes, 1);
    fprintf(fid, '%d\n', nNodes);
    fprintf(fid, '%d %f %f %f\n', nodes.');
    
    fprintf(fid, '$EndNodes\n');
end

function write_elements(fid, elements)
    fprintf(fid, '$Elements\n');
    
    nElems = size(elements, 1);
    fprintf(fid, '%d\n', nElems);
    
    for iElem = 1 : nElems
        fprintf(fid, '%d ', iElem);
        switch elements(iElem, 2)
            % 3-node tria
            case {23, 231}
                fprintf(fid, '%d %d %d %d ', 2, 2, 0, 0);
                fprintf(fid, '%d %d %d\n',elements(iElem, 5:7));
            % 4-node quad
            % TODO: check node ordering
            case {24, 241}
                fprintf(fid, '%d %d %d %d ', 3, 2, 0, 0);
                fprintf(fid, '%d %d %d %d\n',elements(iElem, 5:8));
        end
    end
    
    fprintf(fid, '$EndElements\n');
end

