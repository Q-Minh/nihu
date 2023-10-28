function [tree, m2m] = build_tree(D, depth, nExp, kernel)
%BUILD_TREE build cluster tree for bbFMM method
%   [tree, m2m] = build_tree(D, depth, nExp, kernel) builds the cluster
%   tree and returns the matrix of L2L and M2M translations for the 1D
%   bbFMM method.

tree = struct(...
    'nClus', num2cell(2.^(1:depth)), ...
    'dClus', num2cell(D./2.^(1:depth)), ...
    'level', num2cell(1:depth));

ch0 = chebroots(nExp);

m2m = chebinterp(nExp, [
    lintrans(ch0, [-1 1]', [-1 0]');
    lintrans(ch0, [-1 1]', [0 1]')
    ]);

for level = 2 : depth
    fprintf(1, 'Building level %2d / %2d\n', level, depth);
    
    tree(level).Local = zeros(nExp*tree(level).nClus, 1);
    
    % M2L
    n = tree(level).nClus;
    d = tree(level).dClus;
    
    x0 = lintrans(ch0, [-1 1]', [0 d]');
    
    Ind = zeros(n,3);
    for i = 1 : n/2
        Ind((i-1)*2+(1:2),:) = [-1 3 4; -1 0 4]+2*(i-1);
    end
    
    m = sum(sum(Ind > 0 & Ind <= n));
    
    End = 0;
    I = zeros(m*nExp, nExp);
    J = zeros(m*nExp, nExp);
    K = zeros(m*nExp, nExp);
    for i = 1 : size(Ind,1)
        for j = 1 : size(Ind,2)
            k = Ind(i,j);
            if k < 1 || k > n, continue; end
            I(End+(1:nExp),:) = repmat((i-1)*nExp+(1:nExp)', 1, nExp);
            J(End+(1:nExp),:) = repmat((k-1)*nExp+(1:nExp), nExp, 1);
            K(End+(1:nExp),:) = kernel(x0+i*d, x0+k*d);
            End = End + nExp;
        end
    end
    tree(level).M2L = sparse(I(:), J(:), K(:));
end

end % of function build_tree
