function export_excitation(freqs, exc, fname)

nFreqs = numel(freqs);
nDof = size(exc,1);

fid = fopen(fname, 'w');
fprintf(fid, '%u %u\n', nDof, nFreqs);
fclose(fid);

data = nan(nFreqs, 2*nDof+1);
data(:,1) = freqs;
data(:,2:2:end) = real(exc.');
data(:,3:2:end) = imag(exc.'); %#ok<NASGU>

save(fname, 'data', '-ascii', '-append');
end
