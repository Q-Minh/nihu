function a = import_matrix(fname)

fid = fopen(fname, 'rt');
Us = fscanf(fid, ' (%g,%g)', [2, Inf]);
fclose(fid);
a = complex(Us(1,:), Us(2,:));
n = size(a,2);
n = sqrt(n);
a = reshape(a, n, n).';

end
