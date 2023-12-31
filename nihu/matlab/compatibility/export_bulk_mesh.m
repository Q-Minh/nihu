function export_bulk_mesh(model, fname)
%EXPORT_BULK_MESH Export NiHu mesh to bulk format
%   MODEL = EXPORT_BULK_MESH(FNAME) exports the NiHu mesh into a Nastran
%   Bulk file.
%
% See also: IMPORT_BULK_MESH

%   Copyright 2008-2013 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified: 2017.06.23.

%% Open the file
fid = fopen(fname, 'w');

%% Write HEADER
fprintf(fid, [repmat('$', 1, 72) '\n']);
fprintf(fid, '$\n$ Exported by NiHu (%s.m)\n', mfilename);
fprintf(fid, '$ Created at : %s\n$\n', datestr(now));
fprintf(fid, [repmat('$', 1, 72) '\n']);
fprintf(fid, 'BEGIN BULK\n');

%% Write GRID data
fprintf(fid, 'GRID*   %16d                % 16G% 16G\n*       % 16G\n',...
    model.Nodes.');

%% Write Element data
ElementTypes = {
    'CBAR    ', ShapeSet.LinearLine
    'CTRIA3  ', ShapeSet.LinearTria
    'CQUAD4  ', ShapeSet.LinearQuad
    'CTETRA  ', ShapeSet.LinearTetra
    'CPENTA  ', ShapeSet.LinearPenta
    'CHEXA   ', ShapeSet.LinearHexa
    };
for net = 1 : size(ElementTypes,1)
    LSet = ElementTypes{net, 2};
    elem = model.Elements(:,2) == LSet.Id;
    if any(elem)
        if (LSetId < 100)
            nnod = mod(LSetId, 10);
        else
            nnod = size(LSet.Nodes,1);
        end
        name = ElementTypes{net,1};
        if nnod > 6
            format = [name '%8d       1%8d%8d%8d%8d%8d%8d\n        '...
                repmat('%8d', 1, nnod-6) '\n'];
        else
            format = [name '%8d       1'...
                repmat('%8d', 1, nnod) '\n'];
        end
        fprintf(fid, format, model.Elements(elem,[1 4+(1:nnod)])');
    end
end
%%
fprintf(fid, 'ENDDATA\n');
fclose(fid);
end
