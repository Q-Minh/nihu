function ps = planewave_cyl2d(x, a, k, m)
%PLANEWAVE_CYL2D scattering of an acoustic plane wave from a cylinder
%   ps = PLANEWAVE_CYL2D(X, A, K, M) computes the acoustic wave field
%   scattered from a rigid cylinder of radius A at a wave number K. The
%   scattered field is computed in points described by the rows of matrix X

% Simulate planewave reflection from cylinder
% Copyright 2011-2012 Bence Olteán.

[phi, r] = cart2pol(x(:,1), x(:,2));

mm = 1:m;
ps = zeros(length(r),m+1);

g0 = atan(-besselj(1,k*a)/bessely(1,k*a));	
A0 = -1i*exp(-1i*g0)*sin(g0);
ps(:,1) = repmat(A0,length(r),1).*(besselj(0,k*r)+1i*bessely(0,k*r));

gm = atan((besselj(mm-1,k*a)-besselj(mm+1,k*a))./(bessely(mm+1,k*a)-bessely(mm-1,k*a)));
Am = -2*1i.^(mm+1).*exp(-1i*gm).*sin(gm);
ps(:,2:end) = repmat(Am,length(r),1).*cos(phi*mm).*...
    (besselj(repmat(mm, length(r),1),repmat(k*r, 1, m))+...
    1i*bessely(repmat(mm, length(r), 1),repmat(k*r, 1, m)));

ps = conj(sum(ps,2));

end
