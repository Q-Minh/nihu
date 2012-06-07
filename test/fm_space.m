function [B, C, M] = fm_space(tree, intdata)

depth = length(tree);
B = zeros(length(tree),1);
C = zeros(length(tree),1);
for l = 1 : depth
    B(l) = size(tree(l).D,1) * length(intdata(l).W);
    if isfield(tree(l), 'D2')
        B(l) = B(l) + size(tree(l).D2,1) * length(intdata(l).W);
    end
    C(l) = size(tree(l).coord,1) * length(intdata(l).W);
end

M = sum(C) + C + B;