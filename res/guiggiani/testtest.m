clear;

xi0 = [.66 .66];
[L, dL, ddL] = shapefun(xi0, 242);
m = sqrt(2)/2;
X = [
    1 1 1 m 0 0 0 m m
    0 1 2 2 2 1 0 0 1
    0 0 0 m 1 1 1 m m
    ];
N0 = 1;
N1 = 0;

theta = -2.3365891293288557;

x0 = X * squeeze(L)';
dx0 = X * squeeze(dL);
ddx0 = X * squeeze(ddL);

Avec = dx0(:,1) * cos(theta) + dx0(:,2) * sin(theta);
A = norm(Avec);
J0vec = cross(dx0(:,1), dx0(:,2));

nx = J0vec / norm(J0vec);

Bvec = ddx0(:,1) * cos(theta)^2/2 +...
    ddx0(:,2) * cos(theta) * sin(theta) +...
    ddx0(:,3) * sin(theta)^2 / 2;
J1vec = (cross(ddx0(:,1), dx0(:,2)) + cross(dx0(:,1), ddx0(:,2))) * cos(theta) +...
    (cross(ddx0(:,2), dx0(:,2)) + cross(dx0(:,1), ddx0(:,3))) * sin(theta);

g1vec = Avec/A^2 * (dot(Bvec, J0vec) + dot(Avec, J1vec));

b0vec = -J0vec;
b1vec = 3*g1vec - J1vec;

a0vec = b0vec * N0;
a1vec = b1vec * N0 + b0vec * N1;

Sm2 = -3*dot(Avec, Bvec) / A^5;
Sm3 = 1/A^3;

Fm1 = -(Sm2 * a0vec + Sm3 * a1vec) / (4*pi);
Fm2 = -(Sm3 * a0vec) / (4*pi);

fprintf(1, 'Fm1: %g\n', dot(Fm1, nx));
fprintf(1, 'Fm2: %g\n', dot(Fm2, nx));
