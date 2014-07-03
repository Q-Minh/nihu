function display_block_structure(Ctree, B)

N = length(Ctree(1).ind);
M = zeros(N,N);
m = 0;
for b = 1 : length(B)
    C1 = Ctree(B(b,1));
    C2 = Ctree(B(b,2));
    n = length(C1.ind)*length(C2.ind);
    m = m + n;
    M(C1.ind, C2.ind) = log10(n);
%     M(C1.ind, C2.ind) = b;
end

p = surf(M);
title(sprintf('nz: %d', m));
view(2);
set(p, 'LineStyle', 'none');
set(gca, 'yDir', 'reverse');
axis equal tight;

end
