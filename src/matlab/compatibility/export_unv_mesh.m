function export_unv_mesh(model, fname)
fid = fopen(fname, 'w');
% export nodes
fprintf(fid, '%6d\n', -1);
fprintf(fid, '%6d\n', 2411);
for i = 1 : size(model.Nodes,1)
    fprintf(fid, '%10d%10d%10d%10d\n', model.Nodes(i,1), 1, 1, 1);
    fprintf(fid, '%s\n', ...
        strrep(sprintf('%25.16E%25.16E%25.16E', model.Nodes(i,2:4)), 'E', 'D'));
end
fprintf(fid, '%6d\n', -1);
% export elements
fprintf(fid, '%6d\n', -1);
fprintf(fid, '%6d\n', 2412);
for i = 1 : size(model.Elements,1)
    nNodes = size(ShapeSet.fromId(model.Elements(1,2)).Nodes,1);
    if nNodes == 3
        fedesc = 91;
    elseif nNodes == 4
        fedesc = 94;
    elseif nNodes == 2
        error('Not ready for handling two node elements');
    end
    fprintf(fid, '%10d%10d%10d%10d%10d%10d\n', model.Elements(i,1), fedesc,...
        1, 1, 1, nNodes);
    fprintf(fid, '%10d', model.Elements(i,4+(1:nNodes)));
    fprintf(fid, '\n');
end
fprintf(fid, '%6d\n', -1);
fclose(fid);
end
