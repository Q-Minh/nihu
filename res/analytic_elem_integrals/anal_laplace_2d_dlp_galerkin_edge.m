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
n = [s, -c];
rdn = simple((rvec * n.') / r);

G = simple(-rdn/r/2/pi);

Iinner = simple(int(G, y));
Iinner = simple(subs(Iinner, y, d2) - subs(Iinner, y, eps));
Iinner = simple(limit(Iinner, eps, 0));

Iouter = simple(int(Iinner, x));
Iouter = simple(subs(Iouter, x, d1) - subs(Iouter, x, eps));
Iouter = simple(limit(Iouter, eps, 0));

pretty(Iouter)

%% reconstruction
f11 = d1 * (atan((d2+d1*c)./(d1*s)) - atan(c./s));
f12 = d2 * (atan((d1+d2*c)./(d2*s)) - atan(c./s));
R3 = sqrt(d1^2+d2^2+2*d1*d2*c);
func = -(f11-f12*c + d2*s*log(R3./d2))/(2*pi);

%% difference
difference = Iouter - func;
phi = sym('phi', 'real');
simple(subs(subs(difference, s, sin(phi)), c, cos(phi)))

%% Taylor series
t0 = simple(limit(diff(subs(subs(func, s, sin(phi)), c, cos(phi)), phi, 0), phi, 0))/factorial(0);
t1 = simple(limit(diff(subs(subs(func, s, sin(phi)), c, cos(phi)), phi, 1), phi, 0))/factorial(1);
t2 = simple(limit(diff(subs(subs(func, s, sin(phi)), c, cos(phi)), phi, 2), phi, 0))/factorial(2);
t3 = simple(limit(diff(subs(subs(func, s, sin(phi)), c, cos(phi)), phi, 3), phi, 0))/factorial(3);
