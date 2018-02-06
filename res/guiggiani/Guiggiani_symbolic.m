function [Fm1_ddG, Fm2_ddG, Fm1_dG, Fm2_dG, Fm1_G, Fm2_G] =...
    Guiggiani_symbolic(ndim)

if nargin == 0
    ndim = 3;
end

J0 = sym('J0', [ndim 1]); assume(J0, 'real');
J1 = sym('J1', [ndim 1]); assume(J1, 'real');
N0 = sym('N0', 'real'); assume(N0, 'real');
N1 = sym('N1', 'real'); assume(N1, 'real');
rho = sym('rho', 'real');
nx = sym('nx', [ndim 1]); assume(nx, 'real');
A = sym('A'); assume(A, 'positive');
n = sym('n'); assume(n, 'real');
r1 = sym('r1', [ndim 1]); assume(r1, 'real');
r2 = sym('r2', [ndim 1]); assume(r2, 'real');

d0 = r1/A;
d1 = r2/A - r1*(r1.'*r2)/A^3;

gradr = d0 + rho * d1;
J = J0 + rho * J1;
N = N0 + rho * N1;

assume(r1.'*J0 == 0);

rdny = gradr.'*J;
rdnx = -gradr.'*nx;

rn = rho^n * A^n * (1+n*rho*dot(r1,r2)/A^2);

if (ndim == 3)
    G = subs(rn, n, -1) * N;
    dG = -subs(rn, n, -2) * rdny * N;
    ddG = subs(rn, n, -3) * (3 * rdny * rdnx + (J.' *  nx)) * N;

    G = G * rho;
    dG = dG * rho;
    ddG = ddG * rho;
end


Fm2_ddG = simplify(limit(ddG * rho^2, rho, 0));
Fm1_ddG = simplify(limit(diff(ddG * rho^2, rho, 1), rho, 0));

Fm2_dG = simplify(limit(dG * rho^2, rho, 0));
Fm1_dG = simplify(limit(diff(dG * rho^2, rho, 1), rho, 0));
Fm0_dG = simplify(limit(diff(dG * rho^2, rho, 2)/2, rho, 0));

Fm2_G = simplify(limit(G * rho^2, rho, 0));
Fm1_G = simplify(limit(diff(G * rho^2, rho, 1), rho, 0));
end
