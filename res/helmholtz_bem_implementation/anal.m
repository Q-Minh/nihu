clear;
clc;

syms r R rdnx rdny nx1 nx2 nx3 ny1 ny2 ny3 r1 r2 r3 z theta

G = 1/r;

dG = diff(G, r);
ddG = diff(dG, r);

rdny = (r1*ny1 + r2*ny2 + r3*ny3) / r;
rdnx = -(r1*nx1 + r2*nx2 + r3*nx3) / r;
nxny = nx1*ny1 + nx2*ny2 + nx3*ny3;

% K = G;
% K = dG * rdny;
% K = dG * rdnx;
K = (ddG - dG/r) * rdnx * rdny - dG / r * nxny;

K = subs(K, ny1, 0);
K = subs(K, ny2, 0);
K = subs(K, ny3, 1);
K = subs(K, r, sqrt(R^2 + z^2));
K = subs(K, r1, R * cos(theta));
K = subs(K, r2, R * sin(theta));
K = subs(K, r3, -z);

integrand = K * R;
integral = int(integrand, R);

% pretty(simplify(integral))

syms eps
pretty(simplify(limit(integral - subs(integral, R, eps), eps, 0)))

