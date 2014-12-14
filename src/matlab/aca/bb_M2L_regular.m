function [M2L, trans] = bb_M2L_regular(dim, kernel, nExp)

% generate displacement vectors
trans = transgen(dim, 7);

% preallocate space for the M2L matrices
M2L = zeros(nExp^dim, nExp^dim, size(trans, 1));

% unit source cluster around the origin
I = repmat([-1;1], 1, dim);
sou = lintrans(chebroots(nExp, dim), I, I/2);

for i = 1 : size(trans,1)
    % receiver cluster nodes
    rec = bsxfun(@plus, sou, trans(i,:));
    % compute transfer matrix
    M2L(:,:,i) = kernel(rec, sou);
end

end
