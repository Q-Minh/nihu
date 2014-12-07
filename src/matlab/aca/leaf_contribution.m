function [P2M, L2P] = leaf_contribution(x, y, level, nExp)

Nx = length(x);
Ny = length(y);

End = 0;
I = zeros(Ny,nExp);
J = zeros(Ny,nExp);
K = zeros(Ny,nExp);
father = ceil(y/level.dClus);
for i = 1 : level.nClus
    inter = level.dClus * [i-1, i]';
    index = find(father == i);
    n = length(index);
    I(End+(1:n),:) = repmat(index(:), 1, nExp);
    J(End+(1:n),:) = repmat((i-1)*nExp + (1:nExp), length(index), 1);
    K(End+(1:n),:) = chebinterp(nExp, y(index), inter)';
    End = End + n;
end
P2M = sparse(J(:), I(:), K(:), level.nClus*nExp, Ny);

End = 0;
I = zeros(Nx,nExp);
J = zeros(Nx,nExp);
K = zeros(Nx,nExp);
father = ceil(x/level.dClus);
for i = 1 : level.nClus
    inter = level.dClus * [i-1, i]';
    index = find(father == i);
    n = length(index);
    I(End+(1:n),:) = repmat(index(:), 1, nExp);
    J(End+(1:n),:) = repmat((i-1)*nExp + (1:nExp), length(index), 1);
    K(End+(1:n),:) = chebinterp(nExp, x(index), inter)';
    End = End + n;
end
L2P = sparse(I(:), J(:), K(:), Nx, level.nClus*nExp);

end
