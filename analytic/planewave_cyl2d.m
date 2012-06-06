function [ps] = planewave_cyl2d(r, phi, a, k, m)

% Simulate planewave reflection from cylinder
% Copyright 2011-2012 Bence Olteán.

% TODO: - add help and variable explanation
%       - add variable position for cylinder and plane wave incident angle

mm = 1:m;
ps = zeros(length(r),m+1);

g0 = atan(-besselj(1,k*a)/bessely(1,k*a));	
A0 = -1i*exp(-1i*g0)*sin(g0);
ps(:,1) = repmat(A0,length(r),1).*(besselj(0,k*r)+1i*bessely(0,k*r));

gm = atan((besselj(mm-1,k*a)-besselj(mm+1,k*a))./(bessely(mm+1,k*a)-bessely(mm-1,k*a)));
Am = -2*1i.^(mm+1).*exp(-1i*gm).*sin(gm);
ps(:,2:end) = repmat(Am,length(r),1).*cos(phi*mm).*(besselj(mm,k*r)+1i*bessely(mm,k*r));

ps = sum(ps,2);

end
