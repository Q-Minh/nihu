function P2M = bb_P2M(x, CTree, x0, nExp)
%BB_P2M Point to Multipole sparse matrix

nClusters = length(CTree);
dim = size(x0,2);
nNod = nExp^dim;
nRec = nClusters * nNod;
nSou = size(x,1);

Capacity = 5000;
II = zeros(Capacity,1);
JJ = zeros(Capacity,1);
ZZ = zeros(Capacity,1);

n = 0;
for c = 1 : nClusters
    ch = CTree(c).children;
    if ~isempty(ch)
        continue;
    end
    i = (c-1)*nNod + (1:nNod)';
    j = CTree(c).ind(:);
    
    [J, I] = meshgrid(j,i);
    z = chebinterp(nExp, x(j,:), CTree(c).bb);
    
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
P2M = sparse(II, JJ, ZZ, nRec, nSou);

end

