function fi = bb_far_transfer(tree, sigma, M2M, M2L, L2L)

N = size(M2L,1);
depth = length(tree)-1;

for level = 0 : depth
    iL = level + 1;
    res(iL).multi = zeros(N,size(tree(iL).coord,1));
    res(iL).local = zeros(N,size(tree(iL).coord,1));
end

% compute locals
for level = depth : -1 : 0
    iL = level + 1;
    
    [s, ~, r] = find(tree(iL).interlist);
    [~, ~, ridx] = find(tree(iL).interlistidx);
    for j = 1 : length(s)
        res(iL).local(:,r(j)) = res(iL).local(:,r(j)) + ...
            M2L(:,:,ridx(j)) * res(iL).multi(:,s(j));
    end
end

