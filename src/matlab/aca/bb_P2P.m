function P2P = bb_P2P(B_near, RowTree, ColTree, x, y, kernel)
%BB_P2P Point to Point sparse matrix
%   P2P = bb_P2P(B_NEAR, ROWTREE, COLTREE, X, Y, KERNEL) computes the Point
%   to Point transfer matrix of the black box FMM method.
% B_NEAR  denotes the near field block structure matrix
% ROWTREE  denotes the row cluster tree
% COLTREE  denotes the column cluster tree
% X and Y  denote the receiver and source points, respectively
% KERNEL denotes the matrix kernel
%
% see also: bb_M2L bb_M2M bb_P2M

% Copyright (c) 2014-2014 Peter Fiala

nBlocks = size(B_near,1);

Capacity = 5000;
II = zeros(Capacity,1);
JJ = zeros(Capacity,1);
ZZ = zeros(Capacity,1);

n = 0;
for b = 1 : nBlocks
    block = B_near(b,:);
    
    i = RowTree(block(1)).ind;
    j = ColTree(block(2)).ind;
    [J, I] = meshgrid(j,i);
    z = kernel(x(i,:), y(j,:));
    
    idx = n+(1:numel(J));
    increase = false;
    while max(idx) > Capacity
        increase = true;
        Capacity = Capacity * 2;
    end
    if increase
        II(Capacity) = 0;
        JJ(Capacity) = 0;
        ZZ(Capacity) = 0;
    end

    II(idx) = I(:);
    JJ(idx) = J(:);
    ZZ(idx) = z(:);

    n = n + numel(J);
end
II = II(1:n);
JJ = JJ(1:n);
ZZ = ZZ(1:n);
P2P = sparse(II, JJ, ZZ, size(x,1), size(y,1));

end

