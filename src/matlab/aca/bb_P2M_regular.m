function P2M = bb_P2M_regular(nExp, nClus, fathersou, dsrc)
%BB_P2M_REGULAR Point to Multipole sparse matrix

nPts = size(fathersou,1);
dim = size(dsrc,2);
N = nExp^dim;

II = zeros(nPts*N,1);
JJ = zeros(nPts*N,1);
KK = zeros(nPts*N,1);
n = 0;
for iClus = 1 : nClus
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

end
