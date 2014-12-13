function trans = transgen(dim, n)

idx = (0 : n^dim-1)';
trans = zeros(length(idx),dim);
offset = (n-1)/2;
for d = 1 : dim
    trans(:,d) = mod(idx, n)-offset;
    idx = floor(idx / n);
end

end
