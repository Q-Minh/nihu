clear all;

x = sym('x', 'positive');
y = sym('y', 'positive');
theta = sym('theta', 'real');
c = sym('c', 'real');
s = sym('s', 'real');
d1 = sym('d1', 'positive');
d2 = sym('d2', 'positive');
eps = sym('eps', 'positive');

%% analytical derivation
rvec = [x+y*c, y*s];
r = sqrt(rvec * rvec.');

G = simple(-log(r)/2/pi);

Iinner = simple(int(G, y));
Iinner = simple(subs(Iinner, y, d2) - subs(Iinner, y, eps));
Iinner = simple(limit(Iinner, eps, 0));

Iouter = simple(int(Iinner, x));
Iouter = simple(subs(Iouter, x, d1) - subs(Iouter, x, eps));
Iouter = simple(limit(Iouter, eps, 0));

pretty(Iouter)

%% reconstruction
R3 = sqrt(d1^2 + 2*d1*d2*c + d2^2);
func = (...
    d1*d2 * (3-2*log(R3)) ...
    + s*d1^2*(atan(c/s) - atan((d2 + d1*c)/(d1*s)))...
    + s*d2^2*(atan(c/s) - atan((d1 + d2*c)/(d2*s)))...
    + c*(d1^2*log(d1/R3) + d2^2*log(d2/R3))...
    ) / (4*pi);

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
subs(subs(subs(subs(func, d1, r1), d2, r2), c, cos(phi)), s, sin(phi))