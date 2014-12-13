function [M2L, trans] = bb_M2L_regular(dim, kernel, nExp)

trans = transgen(dim, 7);

I = repmat([-1;1], 1, dim);
sou = lintrans(chebroots(nExp, dim), I, I/2);

M2L = zeros(nExp^dim, nExp^dim, length(trans));
for i = 1 : size(trans,1)
    rec = bsxfun(@plus, sou, trans(i,:));
    M2L(:,:,i) = kernel(rec, sou);
end

end
