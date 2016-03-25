function [U, V] = lowrank_approx(M, siz, eps, R)
%LOWRANK_APPROX low rank approximation of a matrix
%   [U,V] = LOWRANK_APPROX(M, siz, eps, R) returns the low-rank
%   approximation of the matrix M.
%   The low rank approximation is of the form
%   M = U * V'
%   siz denotes the size of the matrix, eps is the required
%   relative error and R is the maximal rank of the approximation.
%
% See also: LOWRANK_APPROX_BLOCK

% Copyright (C) 2014-2015 Peter Fiala

nRows = siz(1); % number of rows
nCols = siz(2); % number of columns

if nargin < 4   % maximal possible rank
    R = min(nRows, nCols);
end

U = zeros(nRows,R); % preallocate space for low rank
V = zeros(nCols,R);
S2 = 0;             % norm of estimator

r = 0;                      % result counter
uncheckedrows = 1:nRows;    % indices of unchecked rows
i = 1;                      % actual investigated row

% estimate magnitude of matrix elements by computing one percent
C = min(ceil(nRows*nCols/100), R);
irand(:,1) = randperm(nRows, C);
irand(:,2) = randperm(nCols, C);
magest = 0;
for k = 1 : C
    magest = magest+abs(M(irand(k,1),irand(k,2)));
end
magest = magest/C;

for k = 1 : R   % max R iterations should be enough
    uncheckedrows = setdiff(uncheckedrows,i);
    
    row = M(i,1:nCols) - U(i,1:r) * V(:,1:r)';
    [gamma, j] = max(abs(row));
    if gamma / magest < 1e-8
        fprintf(1, 'Full zero row: %3d\n', i);
        i = uncheckedrows(1);
        continue;
    end
    
    col = M(1:nRows,j) - U(:,1:r) * V(j,1:r)';
    
    r = r + 1;
    
    U(:, r) = col/col(i);
    V(:, r) = row';
    
    S2 = S2 + norm(U(:,r), 'fro')^2 * norm(V(:,r), 'fro')^2;
    for j = 1 : r-1
        S2 = S2 + 2 * (U(:,r)'*U(:,j)) * (V(:,j)'*V(:,r));
    end
    
    if r > 1
        err = norm(U(:,r), 'fro') * norm(V(:,r), 'fro') / sqrt(S2);
        
        if err < eps
            break;
        end
    end
    
    [~, idx] = max(abs(col(uncheckedrows)));
    i = uncheckedrows(idx);
end % of loop over iterations

U = U(:, 1:r);
V = V(:, 1:r);

end % of function
