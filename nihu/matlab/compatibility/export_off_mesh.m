function export_off_mesh(mesh, fname)
%EXPORT_OFF_MESH Export a NiHu mesh to OFF format

[nodes, elements] = extract_core_mesh(mesh);
% delete trailing empty columns
elements = elements(:,1:find(any(elements,1), 1, 'last'));

elements(elements(:,1) == 32404, 1) = 4;
elements(elements(:,1) == 32303, 1) = 3;
elements(elements(:,1) == 21202, 1) = 2;

elements(elements(:,1) == 21203, 1) = 3;

fid = fopen(fname, 'wt');
if fid == -1
    error('NiHu:file_open', 'Cannot open file %s', fname);
end

fprintf(fid, 'OFF\n');
fprintf(fid, '%d %d 0\n', size(nodes,1), size(elements,1));
fprintf(fid, '%g %g %g\n', nodes.');
for e = 1 : size(elements,1)
    n = elements(e,1);
    fprintf(fid, [repmat('%u ', 1, n+1) '\n'], elements(e,1:n+1));
end

fclose(fid);
end
