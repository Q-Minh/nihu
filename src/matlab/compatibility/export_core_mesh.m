function export_core_mesh(mesh, fname)
%EXPORT_CORE_MESH export a NiHu mesh to simple text format
%   EXPORT_CORE_MESH(MESH, FNAME) exports the NiHu mesh MESH into the file
%   FNAME.

% see also EXPORT_EXCITATION EXTRACT_CORE_MESH

[nodes, elements] = extract_core_mesh(mesh);
% delete trailing empty columns
elements = elements(:,1:find(any(elements,1), 1, 'last'));

fid = fopen(fname, 'wt');
if fid == -1
    error('NiHu:file_open', 'Cannot open file %s', fname);
end

fprintf(fid, '%d\n', size(nodes,1));
fprintf(fid, '%g %g %g\n', nodes.');

fprintf(fid, '%d\n', size(elements,1));
fprintf(fid, [repmat('%u ', 1, size(elements,2)) '\n'], elements.');

fclose(fid);
end