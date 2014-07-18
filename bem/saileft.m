function M = saileft(A, S)
%SAILEFT Sparse approximate inverse
%   M = SAILEFT(A, S) computes the sparse approximate inverse of A.
% Parameters:
%   A : N x N sparse matrix to be inverted
%   S : Sparsity pattern of the inverse M. Two methods are supported.
%       If size(S,1) = N, then M will have the same sparsity pattern as S.
%       If size(S,1) = B < N, then M will be block diagonal with B full
%         blocks. The element indices of the b-th block are the nonzero
%         entries of S(b,:).
%       The default value of S is  A ~= 0
%
% See also: saistruct, sairight

% Peter Fiala
% 2009

%% Parameter check
if nargin == 1
    S = A ~= 0;
end
N = size(A,1);
B = size(S,1);

%% Computing inverse
if B == N % Method 1
    % allocate space for sparese matrix
    nnonzero = length(find(S));
    rowind = zeros(nnonzero,1);
    colind = zeros(nnonzero,1);
    values = zeros(nnonzero,1);
    % starting indices of each new row
    limits = [0; cumsum(sum(S ~= 0, 2))];
    for n = 1 : N
        ind = limits(n)+1:limits(n+1);
        j = find(S(n,:));
        a = A(j,:);
        m = a(:,n)' / (a * a'); % LS solution;
        progbar(1,N,n);
        rowind(ind) = n;
        colind(ind) = j;
        values(ind) = m(:);
        progbar(limits(2),limits(end),limits(n+1));
    end
    M = sparse(rowind, colind, values);
else % Method 2
    blocksize = sum(S~=0, 2);
    nnonzero = sum(blocksize.^2);
    rowind = zeros(nnonzero,1);
    colind = zeros(nnonzero,1);
    values = zeros(nnonzero,1);
    ind = 0;
    for b = 1 : B
        j = S(b,:);     % block indices
        j = j(j ~= 0);  % without zeros
        nj = length(j); % block size
        % compute full inverse of block
        a = A(j,:);
        m = a(:,j)' / (a * a');  % LS solution
        % place the result into the sparse matrix
        [ci, ri] = meshgrid(j,j);
        ind = max(ind) + (1:nj^2);
        rowind(ind) = ri(:);
        colind(ind) = ci(:);
        values(ind) = m(:);
        progbar(1,B,b);
    end
    M = sparse(rowind, colind, values);
end
end
