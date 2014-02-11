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
n = [s, -c];
rdn = simple((rvec * n.') / r);

G = simple(rdn/r);

Iinner = simple(int(G, y));
Iinner = simple(subs(Iinner, y, R2) - subs(Iinner, y, eps));
Iinner = simple(limit(Iinner, eps, 0));

Iouter = simple(int(Iinner, x));
Iouter = simple(subs(Iouter, x, R1) - subs(Iouter, x, eps));
Iouter = simple(limit(Iouter, eps, 0));

pretty(Iouter)

%% reconstruction
f11 = R1 * (atan((R2+R1*c)./(R1*s)) - atan(c./s));
f12 = R2 * (atan((R1+R2*c)./(R2*s)) - atan(c./s));
R3 = sqrt(R1^2+R2^2+2*R1*R2*c);
func = f11-f12*c + R2*s*log(R3./R2);

%% difference
difference = Iouter - func;
phi = sym('phi', 'real');
simple(subs(subs(difference, s, sin(phi)), c, cos(phi)))

%% Taylor series
t0 = simple(limit(diff(subs(subs(func, s, sin(phi)), c, cos(phi)), phi, 0), phi, 0))/factorial(0);
t1 = simple(limit(diff(subs(subs(func, s, sin(phi)), c, cos(phi)), phi, 1), phi, 0))/factorial(1);
t2 = simple(limit(diff(subs(subs(func, s, sin(phi)), c, cos(phi)), phi, 2), phi, 0))/factorial(2);
t3 = simple(limit(diff(subs(subs(func, s, sin(phi)), c, cos(phi)), phi, 3), phi, 0))/factorial(3);
