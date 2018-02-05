clear all;

nu = sym('nu', 'real');
mu = sym('mu', 'real');

rho = sym('rho', 'positive');
eps = sym('eps', 'positive');
R = sym('R', 'positive');
z = sym('z', 'real');
theta = sym('theta', 'real');
T = sym('T', 'real');

rvec = [rho*cos(theta); rho*sin(theta); z];
r = norm(rvec);
rdy = rvec/r;

A = 1/(16*pi*mu*(1-nu)*r);
U = A * ((3-4*nu)*eye(3) + (rdy*rdy.'));

U0 = simple(subs(U, z, 0));
Iinner = simple(int(U0*rho, rho));
Iinner = simple(subs(Iinner, rho, R/cos(theta))-subs(Iinner, rho, eps));
Iinner = limit(Iinner, eps, 0);
I = simple(int(Iinner, theta));
I = subs(I, theta, T);
