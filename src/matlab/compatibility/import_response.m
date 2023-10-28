function [p, k] = import_response(fname)
%EXPORT_EXCITATION export a NiHu acoustic BEM excitation into a text file
%   EXPORT_EXCITATION(Q, K, FNAME) exports the NiHu excitation Q at wave
%   number K to the text file FNAME.

fid = fopen(fname, 'rt');
if fid == -1
    error('NiHu:file_open', 'Cannot open file %s', fname);
end

k = fscanf(fid, '%g\n', 1);
s = fscanf(fid, '%u\n', 1);
data = fscanf(fid, '%g %g\n', [2, s]);
p = complex(data(1,:), data(2,:)).';
fclose(fid);
end