function export_core_mesh(mesh, fname)

[nodes, elements] = extract_core_mesh(mesh);
elements = elements(:,1:find(any(elements,1), 1, 'last'));

fid = fopen(fname, 'wt');

fprintf(fid, '%d\n', size(nodes,1));
fprintf(fid, '%g %g %g\n', nodes.');

fprintf(fid, '%d\n', size(elements,1));
fprintf(fid, [repmat('%u ', 1, size(elements,2)) '\n'], elements.');

fclose(fid);
