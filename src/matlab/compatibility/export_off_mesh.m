function export_off_mesh(mesh, fname)
%EXPORT_OFF_MESH Export a NiHu mesh to OFF format

[nodes, elements] = extract_core_mesh(mesh);
% delete trailing empty columns
elements = elements(:,1:find(any(elements,1), 1, 'last'));

elements(elements == 32404) = 4;

fid = fopen(fname, 'wt');
if fid == -1
    error('NiHu:file_open', 'Cannot open file %s', fname);
end

fprintf(fid, 'OFF\n');
fprintf(fid, '%d %d 0\n', size(nodes,1), size(elements,1));
fprintf(fid, '%g %g %g\n', nodes.');
fprintf(fid, [repmat('%u ', 1, size(elements,2)) '\n'], elements.');

fclose(fid);
end
