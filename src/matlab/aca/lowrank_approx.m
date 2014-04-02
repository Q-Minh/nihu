function [U, V] = lowrank_approx(M, I, J, R)

m = length(I);
n = length(J);

if m == 1
    U = 1;
    V = M(I,J);
    return;
end

if n == 1
    U = M(I,J);
    V = 1;
    return;
end

R = min(R, min(m,n));

U = zeros(m,R);
V = zeros(R,n);

% if m <= n
for r = 1 : R
    i = r;
    row = M(I(i),J);
    row = row - U(i,1:r-1) * V(1:r-1,:);
    [~,j] = max(abs(row));
    col = M(I,J(j));
    col = col - U(:,1:r-1) * V(1:r-1,j);
    
    U(:, i) = col/col(i);
    V(i, :) = row;
end
% else
%     for r = 1 : R
%         j = r;
%         col = M(I,J(j));
%         col = col - U(:,1:r-1) * V(1:r-1,j);
%         [~,i] = max(abs(col));
%         row = M(I(i),J);
%         row = row - U(i,1:r-1) * V(1:r-1,:);
%
%         U(:, i) = col/col(i);
%         V(i, :) = row;
%     end
% end

end
