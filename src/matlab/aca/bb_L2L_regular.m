function [L2L, trans] = bb_L2L_regular(dim, nExp)

trans = transgen(dim, 2);

I = repmat([-1;1],1,dim);
xs = lintrans(chebroots(nExp, dim), I, repmat([-.5;.5], 1, dim));

L2L = zeros(nExp^dim, nExp^dim, 2^dim);
for i = 1 : size(trans,1)
    rec = bsxfun(@plus, xs, trans(i,:));
    L2L(:,:,i) = chebinterp(nExp, rec, I)';
end

end
