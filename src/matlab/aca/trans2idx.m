function idx = trans2idx(trans, n)
dim = size(trans,2);
trans = round(trans + (n-1)/2);
powers = n.^(0:dim-1)';
idx = trans * powers + 1;
end
