function M2L = bb_M2L(B_far, x0, y0, rExp, cExp, kernel)
%BB_M2L Multipole to Local sparse matrix
%   M2L = bb_M2L(B_FAR, X0, Y0, REXP, CEXP, KERNEL) computes the multipole
%   to local transfer in a black box FMM computation.
% B_FAR is the far field block structure
% X0 and Y0  denote the receiver and source clusters' Chebyshev nodes,
%   respectively.
% REXP and CEXP denote the Chebyshev expansion of the receiver and source
%   cluster trees
% KERNEL denotes the kernel function
%
% see also: bb_M2M bb_P2M bb_P2P

% copyright (c) 2014-2014 Peter Fiala

nBlocks = size(B_far,1);
dim = size(x0,2);

nCol = cExp^dim;
nRow = rExp^dim;
nEl = nRow*nCol;
N = nBlocks*nEl;
II = zeros(N,1);
JJ = zeros(N,1);
ZZ = zeros(N,1);

n = 0;
for b = nBlocks : -1 : 1
    block = B_far(b,:);
    
    i = (block(1)-1)*nRow + (1:nRow);
    j = (block(2)-1)*nCol + (1:nCol);
    [J, I] = meshgrid(j,i);
    z = kernel(x0(i,:),y0(j,:));
    
    idx = n+(1:nEl);
    II(idx) = I(:);
    JJ(idx) = J(:);
    ZZ(idx) = z(:);
    
    n = n + nEl;
end
M2L = sparse(II, JJ, ZZ, size(x0,1), size(y0,1));

end

