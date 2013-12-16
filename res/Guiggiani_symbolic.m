clear;

d0 = sym('d0', [3,1]); d0 = sym(d0, 'real');
d1 = sym('d1', [3 1]); d1 = sym(d1, 'real');
Jac0 = sym('Jac0', [3 1]); Jac0 = sym(Jac0, 'real');
Jac1 = sym('Jac1', [3 1]); Jac1 = sym(Jac1, 'real');
N0 = sym('N0', 'real'); N0 = sym(N0, 'real');
N1 = sym('N1', 'real'); N1 = sym(N1, 'real');
rho = sym('rho', 'positive');

gradr = d0 + rho * d1;
Jac = Jac0 + rho * Jac1;
N = N0 + rho * N1;

syms n A
Avec = sym('Avec', [3 1]); Avec = sym(Avec, 'real');
Bvec = sym('Bvec', [3 1]); Bvec = sym(Bvec, 'real');

rn = rho^n * A^n * (1+n*rho*dot(Avec,Bvec)/A^2);

nx = sym('nx', [3 1]); nx = sym(nx, 'real');

ddG = (subs(rn, n, -3) * (3 * dot(gradr,Jac) * dot(-gradr, nx) + dot(Jac, nx)) * N * rho);

Fm2 = limit(ddG * rho^2, rho, 0);
Fm1 = limit(diff(ddG * rho^2, rho, 1), rho, 0);
