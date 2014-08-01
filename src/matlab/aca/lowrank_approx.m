function [U, V] = lowrank_approx(M, siz, eps, R)
%LOWRANK_APPROX return low rank approximation of a matrix
%   [U,V] = LOWRANK_APPROX(M, eps) returns the low-rank approximation
%   of the matrix M. The low rank approximation is of the form
%   M = U * V'
%   Argument R denotes the maximal rank of the approximation.
%
% See also: LOWRANK_APPROX_BLOCK
%
% Copyright (C) 2014 Peter Fiala

nRows = siz(1); % number of rows
nCols = siz(2); % number of columns

fprintf('block size: %3dx%3d\n', nRows, nCols);

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
    
    fprintf(1, 'skeleton (%2d,%2d), r=%2d', i, j, r);
    
    S2 = S2 + norm(U(:,r), 'fro')^2 * norm(V(:,r), 'fro')^2;
    for j = 1 : r-1
        S2 = S2 + 2 * (U(:,r)'*U(:,j)) * (V(:,j)'*V(:,r));
    end
    
    if r > 1
        err = norm(U(:,r), 'fro') * norm(V(:,r), 'fro') / sqrt(S2);
        fprintf(1, ', log10 err=%.2f\n', log10(err));
        
        if err < eps
            break;
        end
    else
        fprintf(1, '\n');
    end
    
    [~, idx] = max(abs(col(uncheckedrows)));
    i = uncheckedrows(idx);
end

U = U(:, 1:r);
V = V(:, 1:r);

fprintf(1, 'Compression: %g\n\n', ((nRows+nCols)*r)/(nRows*nCols));

end
