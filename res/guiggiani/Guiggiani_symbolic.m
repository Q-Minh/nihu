function [Fm1_ddG, Fm2_ddG, Fm1_dG, Fm2_dG, Fm1_G, Fm2_G] =...
    Guiggiani_symbolic(ndim)

if nargin == 0
    ndim = 3;
end

% d0 = sym('d0', [ndim,1]); assume(d0, 'real');
% d1 = sym('d1', [ndim 1]); assume(d1, 'real');
Jac0 = sym('Jac0', [ndim 1]); assume(Jac0, 'real');
Jac1 = sym('Jac1', [ndim 1]); assume(Jac1, 'real');
N0 = sym('N0', 'real'); assume(N0, 'real');
N1 = sym('N1', 'real'); assume(N1, 'real');
rho = sym('rho', 'positive');
nx = sym('nx', [ndim 1]); assume(nx, 'real');
A = sym('A'); assume(A, 'positive');
syms n
Avec = sym('Avec', [ndim 1]); assume(Avec, 'real');
Bvec = sym('Bvec', [ndim 1]); assume(Bvec, 'real');

d0 = Avec/A;
d1 = Bvec/A - Avec*dot(Avec, Bvec)/A^3;

gradr = d0 + rho * d1;
Jac = Jac0 + rho * Jac1;
N = N0 + rho * N1;

rn = rho^n * A^n * (1+n*rho*dot(Avec,Bvec)/A^2);

if (ndim == 3)
    G = subs(rn, n, -1) * N;
    dG = -subs(rn, n, -2) * dot(gradr,Jac) * N;
    ddG = subs(rn, n, -3) * (3 * dot(gradr,Jac) * dot(-gradr, nx) + dot(Jac, nx)) * N;

    G = G * rho;
    dG = dG * rho;
    ddG = ddG * rho;
end

Fm2_ddG = simplify(limit(ddG * rho^2, rho, 0));
Fm1_ddG = simplify(limit(diff(ddG * rho^2, rho, 1), rho, 0));

Fm2_dG = simplify(limit(dG * rho^2, rho, 0));
Fm1_dG = simplify(limit(diff(dG * rho^2, rho, 1), rho, 0));

Fm2_G = simplify(limit(G * rho^2, rho, 0));
Fm1_G = simplify(limit(diff(G * rho^2, rho, 1), rho, 0));
end
