clear all;

%% parameters
rho = sym('rho', 'positive');
R = sym('R', 'positive');
theta = sym('theta', 'real');
nu = sym('nu', 'positive');
mu = sym('mu', 'positive');
eps = sym('eps', 'positive');

%% computed bricks
rvec = [rho*cos(theta); rho*sin(theta); 0];
r = simple(norm(rvec));
rdy = simple(rvec/r);
n = [0 0 1].';
rdn = rdy.'*n;

%%
A = 1/(16*pi*mu*(1-nu)*r);
U = simple(A * ((3-4*nu)*eye(3) + (rdy*rdy.')));

%%
Iinner = int(U*rho, rho);
Iinner = subs(Iinner, rho, R/cos(theta)) - subs(Iinner, rho, eps);

Iouter = int(Iinner, theta);
I = simple(limit(Iouter, eps, 0));

%% check
q = 2*atanh(tan(theta/2));
I2 = R * [
    sin(theta)+(3-4*nu)*q, -(cos(theta)+1), 0
    -(cos(theta)+1), 4*(1-nu)*q-sin(theta), 0
    0, 0, (3-4*nu)*q
    ] / (16*pi*mu*(1-nu));

%% the traction kernel
B = 1 / (8*pi*(1-nu)*r^2);
T = simple(-rdn * ((1-2*nu)*eye(3) + 3*(rdy*rdy.')) + (1-2*nu)*(rdy*n.'-n*rdy.'))*B;

%%
Iinner = int(T*rho, rho);
Iupper = subs(Iinner, rho, R/cos(theta));
Ilower = subs(Iinner, rho, eps);

Iu = int(Iupper, theta);
Il = limit(int(Ilower, theta, 0, 2*pi), eps, 0);

Iouter = int(Iinner, theta);
% I = simple(limit(Iouter, eps, 0));

