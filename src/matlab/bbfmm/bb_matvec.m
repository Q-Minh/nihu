function x = bb_matvec(P2P, P2M, M2L, M2M, L2L, L2P, y)
%BB_MATVEC Matrix-vector multiplication with bb FMM
%   X = BB_MATVEC(P2P, P2M, M2L, M2M, L2L, L2P, Y)  performs matrix-vector
%   multiplictaion.

% direct response at low level
x = P2P * y;

%upward and transfer pass
up = P2M * y;
dn = zeros(size(M2L,1),1);
while any(up)
    dn = dn + M2L * up;
    up = M2M * up;
end

% downward pass
while any(dn)
    x = x + L2P * dn;
    dn = L2L * dn;
end

end
