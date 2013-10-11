function export_excitation(q, k, fname)

fid = fopen(fname, 'wt');

fprintf(fid, '%g\n', k);
fprintf(fid, '%u\n', size(q,1));
fprintf(fid, '%g %g\n', [real(q) imag(q)].');

fclose(fid);
end