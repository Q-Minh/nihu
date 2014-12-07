function x = chebroots(n, d)
%CHEBROOTS roots of Chebyshev polynomials
%  X = CHEBROOTS(N) computes the roots of the Chebyshev polynomial Tn(x)
%  X = CHEBROOTS(N, D) returns the roots in D dimensions (1, 2 or 3)
%
% Example:
%   chebroots(5)
%   x = chebroots(5,2);
%   figure; plot(x(:,1), x(:,2), '.');
%   x = chebroots(5,3);
%   figure; plot3(x(:,1), x(:,2), x(:,3), '.');

if nargin < 2
    d = 1;
end

x = cos(pi/2*(2*(n:-1:1)'-1)/n);

if d == 2
    [x1, x2] = meshgrid(x,x);
    x = [x1(:) x2(:)];
elseif d == 3
    [x1, x2, x3] = ndgrid(x,x,x);
    x = [x1(:) x2(:) x3(:)];
elseif d == 4
    [x1, x2, x3, x4] = ndgrid(x,x,x,x);
    x = [x1(:) x2(:) x3(:) x4(:)];
end

end

