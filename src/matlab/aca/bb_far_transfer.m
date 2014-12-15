function [fi, times] = bb_far_transfer(tree, sigma, P2M, M2M, M2L, L2L, L2P, alpha)

N = size(M2L,1);
depth = length(tree)-1;

res = struct(...
    'multi', cell(depth+1, 1), ...
    'local', cell(depth+1, 1));

times = struct(...
    'M2M', cell(depth+1, 1), ...
    'M2L', cell(depth+1, 1), ...
    'L2L', cell(depth+1, 1));

for level = 0 : depth
    iL = level + 1;
    res(iL).multi = zeros(N,size(tree(iL).coord,1));
    res(iL).local = zeros(N,size(tree(iL).coord,1));
end

% insert sigma into leaf level
res(depth+1).multi(:) = P2M*sigma;

% upward pass (M2M)
for level = depth : -1 : 1
    iL = level + 1;
    
    t0 = tic();
    r = tree(iL).father;
    ridx = tree(iL).fatheridx;
    for s = 1 : size(tree(iL).coord,1)
        res(iL-1).multi(:,r(s)) = res(iL-1).multi(:,r(s)) + ...
            M2M(:,:,ridx(s)) * res(iL).multi(:,s);
    end
    times(iL).M2M = toc(t0);
end

% transfer pass (M2L)
for level = 2 : depth
    iL = level + 1;
    
    t0 = tic();
    [s, ~, r] = find(tree(iL).interlist);
    [~, ~, ridx] = find(tree(iL).interlistidx);
    loc = compute_M2L(int32(s)-1, int32(r)-1, int32(ridx)-1, M2L, res(iL).multi);
    res(iL).local = res(iL).local + tree(iL).diameter^alpha * loc;
    %    for j = 1 : length(s)
    %        res(iL).local(:,r(j)) = res(iL).local(:,r(j)) + ...
    %            M2L(:,:,ridx(j)) * tree(iL).diameter^alpha * res(iL).multi(:,s(j));
    %    end
    times(iL).M2L = toc(t0);
end


% downward pass (L2L)
for level = 2 : depth-1
    iL = level + 1;
    
    t0 = tic();
    [s, ~, r] = find(tree(iL).children);
    [~, ~, ridx] = find(tree(iL).childrenidx);
    for j = 1 : length(s)
        res(iL+1).local(:,r(j)) = res(iL+1).local(:,r(j)) + ...
            L2L(:,:,ridx(j)) * res(iL).local(:,s(j));
    end
    times(iL).L2L = toc(t0);
end

fi = L2P * res(end).local(:);

end
