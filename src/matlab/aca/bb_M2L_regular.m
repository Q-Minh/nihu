function [M2L, trans] = bb_M2L_regular(dim, kernel, nExp)

idx = (0 : 7^dim-1)';
trans = zeros(length(idx),dim);
for d = 1 : dim
    trans(:,d) = mod(idx, 7)-3;
    idx = floor(idx / 7);
end

sou = lintrans(chebroots(nExp, dim), repmat([-1;1], 1, dim), repmat([0;1], 1, dim));

M2L = zeros(nExp^dim, nExp^dim, length(trans));
for i = 1 : size(trans,1)
    rec = bsxfun(@plus, sou, trans(i,:));
    M2L(:,:,i) = kernel(rec, sou);
end

end
