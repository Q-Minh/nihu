function [M2M, trans] = bb_M2M_regular(dim, nExp)
%BB_M2M_REGULAR Multipole to Multipole sparse matrix

% Copyright (C) 2014 Peter Fiala

trans = transgen(dim, 2);

I = repmat([-1;1],1,dim);
xs = lintrans(chebroots(nExp, dim), repmat([-1;1], 1, dim), repmat([-.5;.5], 1, dim));

M2M = zeros(nExp^dim, nExp^dim, 2^dim);
for i = 1 : size(trans,1)
    sou = bsxfun(@minus, xs, trans(i,:));
    M2M(:,:,i) = chebinterp(nExp, sou, I);
end

end
