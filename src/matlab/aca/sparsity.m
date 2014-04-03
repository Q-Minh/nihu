function [S, B_far] = sparsity(Ctree, Btree)

Capacity = 5000;
S = zeros(Capacity,2);

N = 0;
keep = true(length(Btree),1);
for b = 1 : length(Btree)
    bt = Btree(b,:);
    i = Ctree(bt(1)).ind;
    j = Ctree(bt(2)).ind;
    if length(i) < 2 || length(j) < 2
        keep(b) = false;
        [I, J] = meshgrid(i, j);
        idx = N + (1:numel(I));
        while max(idx) > Capacity
            Capacity = 2 * Capacity;
            S(Capacity,2) = 0;
        end
        S(idx,1) = I(:);
        S(idx,2) = J(:);
        N = N + numel(I);
    end
end
S = S(1:max(idx),:);

B_far = Btree(keep,:);

end
