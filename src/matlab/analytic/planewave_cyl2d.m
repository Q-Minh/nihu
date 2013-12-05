function ps = planewave_cyl2d(x, a, k, m)
%PLANEWAVE_CYL2D scattering of an acoustic plane wave from a cylinder
%   ps = PLANEWAVE_CYL2D(X, A, K, M) computes the acoustic wave field
%   scattered from a rigid cylinder of radius A at a wave number K. The
%   scattered field is computed in points described by the rows of matrix X

% Simulate planewave reflection from cylinder
% Copyright 2011-2012 Bence Olteán.

[phi, r] = cart2pol(x(:,1), x(:,2));

g = [
    atan(-besselj(1,k*a)/bessely(1,k*a)), ...
    atan((besselj((1:m)-1,k*a)-besselj((1:m)+1,k*a))./(bessely((1:m)+1,k*a)-bessely((1:m)-1,k*a)))
    ];
A = -1i*[1, 2*1i.^(1:m)].*exp(-1i*g).*sin(g);

ps = repmat(A,length(r),1).*cos(phi*(0:m)).*...
    besselh(repmat(0:m, length(r),1),...
    ones(length(r), m+1), repmat(k*r, 1, m+1));

ps = conj(sum(ps,2));

end
