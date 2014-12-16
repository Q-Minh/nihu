function P2M = bb_P2M_regular(nExp, nClus, fathersou, dsrc)
%BB_P2M_REGULAR Point to Multipole sparse matrix

% Copyright (C) 2014 Peter Fiala

nPts = length(fathersou);
dim = size(dsrc,2);
N = nExp^dim;

Z = chebinterp(nExp, dsrc);
II = bsxfun(@plus, (1:N)', N*(fathersou'-1));
JJ = repmat(1:nPts, N, 1);
P2M = sparse(II(:), JJ(:), Z(:), N*nClus, nPts);

end
