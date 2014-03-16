function export_excitation(q, k, fname)
%EXPORT_EXCITATION export a NiHu acoustic BEM excitation into a text file
%   EXPORT_EXCITATION(Q, K, FNAME) exports the NiHu excitation Q at wave
%   number K to the text file FNAME.

% see also EXPORT_EXCITATION EXTRACT_CORE_MESH

fid = fopen(fname, 'wt');
if fid == -1
    error('NiHu:file_open', 'Cannot open file %s', fname);
end

fprintf(fid, '%.8x\n', k);
fprintf(fid, '%u\n', size(q,1));
fprintf(fid, '%.8x %.8x\n', [real(q) imag(q)].');

fclose(fid);
end