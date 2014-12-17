function [fi, times] = bb_far_transfer(tree, sigma, P2M, M2M, M2L, L2L, L2P, alpha)

N = size(M2M,1);
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
    res(iL-1).multi = compute_M2M(res(iL).multi,...
        int32(tree(iL).father)-1,...
        int32(tree(iL).fatheridx)-1,...
        M2M,...
        int32(size(tree(iL-1).coord,1)));
    times(iL).M2M = toc(t0);
end

% transfer pass (M2L)
for level = 2 : depth
    iL = level + 1;
    
    t0 = tic();
    [s, ~, r] = find(tree(iL).interlist);
    [~, ~, ridx] = find(tree(iL).interlistidx);
if isstruct(M2L)
    loc = compute_CM2L(int32(s)-1, int32(r)-1, int32(ridx)-1,...
        M2L.U, M2L.V, int32(M2L.r),...
        res(iL).multi);
else
    loc = compute_M2L(int32(s)-1, int32(r)-1, int32(ridx)-1, M2L, res(iL).multi);
end
    res(iL).local = res(iL).local + tree(iL).diameter^alpha * loc;
    times(iL).M2L = toc(t0);
end


% downward pass (L2L)
for level = 2 : depth-1
    iL = level + 1;
    
    t0 = tic();
    [s, ~, r] = find(tree(iL).children);
    [~, ~, ridx] = find(tree(iL).childrenidx);
    res(iL+1).local = res(iL+1).local + compute_L2L(res(iL).local,...
        int32(s)-1, int32(r)-1, int32(ridx)-1,...
        L2L,...
        int32(size(tree(iL+1).coord,1)));
    times(iL).L2L = toc(t0);
end

fi = L2P * res(end).local(:);

end
