function P2M = bb_P2M_regular(nExp, nClus, fathersou, dsrc)
%BB_P2M_REGULAR Point to Multipole sparse matrix

% Copyright (C) 2014 Peter Fiala

nPts = length(fathersou);
dim = size(dsrc,2);
N = nExp^dim;

<<<<<<< HEAD
II = zeros(nPts*N,1);
JJ = zeros(nPts*N,1);
KK = zeros(nPts*N,1);
n = 0;
for iClus = 1 : nClus
    progbar(1, nClus, iClus);
    idx = (iClus-1)*N + (1:N);
    s = find(fathersou == iClus);
    xs = dsrc(s,:);
    Z = chebinterp(nExp, xs);
    [J, I] = meshgrid(s, idx);
    y = n+(1:numel(J));
    II(y) = I(:);
    JJ(y) = J(:);
    KK(y) = Z(:);
    n = y(end);
end

P2M = sparse(II, JJ, KK, N*nClus, nPts);
=======
Z = chebinterp(nExp, dsrc);
II = bsxfun(@plus, (1:N)', N*(fathersou'-1));
JJ = repmat(1:nPts, N, 1);
P2M = sparse(II(:), JJ(:), Z(:), N*nClus, nPts);
>>>>>>> b95d07095cd9dfb87b6356d922e81e6e58e3d4f8

end
