function [tr, rr, rs, ns] = simpletree(tree, rcoord, rfather, scoord, sfather, snorm, symm)

if nargin < 7
    symm = 0;
end

c = cell(length(tree),1);
tr = struct('father', c, 'interlist', c, 'r', c);

for l = 1 : length(tree)
    tr(l).father = tree(l).father;
    tr(l).interlist = tree(l).interlist.';
    if symm
        tr(l).iminterlist = tree(l).iminterlist.';
    end
    if l > 1
        tr(l).r = (tree(l-1).coord(tree(l).father,:) - tree(l).coord).';
    end
end

rs = (tree(end).coord(sfather,:) - scoord).';
rr = (tree(end).coord(rfather,:) - rcoord).';
ns = snorm.';