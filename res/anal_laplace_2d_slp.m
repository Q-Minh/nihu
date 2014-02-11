clear;

x = sym('x', 'positive');
y = sym('y', 'positive');
theta = sym('theta', 'real');
c = sym('c', 'real');
s = sym('s', 'real');
R1 = sym('R1', 'positive');
R2 = sym('R2', 'positive');
eps = sym('eps', 'positive');

%% analytical derivation
rvec = [x+y*c, y*s];
r = sqrt(rvec * rvec.');

G = simple(log(r));

Iinner = simple(int(G, y));
Iinner = simple(subs(Iinner, y, R2) - subs(Iinner, y, eps));
Iinner = simple(limit(Iinner, eps, 0));

Iouter = simple(int(Iinner, x));
Iouter = simple(subs(Iouter, x, R1) - subs(Iouter, x, eps));
Iouter = simple(limit(Iouter, eps, 0));

pretty(Iouter)

%% reconstruction
R3 = sqrt(R1^2 + 2*R1*R2*c + R2^2);
func = -1/2*(...
    R1*R2 * (3-2*log(R3)) ...
    + s*R1^2*(atan(c/s) - atan((R2 + R1*c)/(R1*s)))...
    + s*R2^2*(atan(c/s) - atan((R1 + R2*c)/(R2*s)))...
    + c*(R1^2*log(R1/R3) + R2^2*log(R2/R3))...
    );

%% difference
difference = Iouter - func;
phi = sym('phi', 'real');
simple(subs(subs(difference, s, sin(phi)), c, cos(phi)))

%% Taylor series
t0 = simple(limit(diff(subs(subs(func, s, sin(phi)), c, cos(phi)), phi, 0), phi, 0))/factorial(0);
t1 = simple(limit(diff(subs(subs(func, s, sin(phi)), c, cos(phi)), phi, 1), phi, 0))/factorial(1);
t2 = simple(limit(diff(subs(subs(func, s, sin(phi)), c, cos(phi)), phi, 2), phi, 0))/factorial(2);
t3 = simple(limit(diff(subs(subs(func, s, sin(phi)), c, cos(phi)), phi, 3), phi, 0))/factorial(3);

%%
phi = pi/2;
r1 = 1;
r2 = 1;
subs(subs(subs(subs(func, R1, r1), R2, r2), c, cos(phi)), s, sin(phi))