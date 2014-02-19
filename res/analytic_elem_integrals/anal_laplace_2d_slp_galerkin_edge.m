clear all;

x = sym('x', 'positive');   % coordinate on x elem
y = sym('y', 'positive');   % coordinate on y elem
c = sym('c', 'real');       % cosine of angle
s = sym('s', 'real');       % sine of angle
d1 = sym('d1', 'positive'); % length of x elem
d2 = sym('d2', 'positive'); % length of y elem
eps = sym('eps', 'positive');   % limiting variable

%% analytical derivation
rvec = [x+y*c, y*s];        % distance vector
r = sqrt(rvec * rvec.');    % scalar distance

G = simple(-log(r)/2/pi);   % Green's function

% primitive function
Iinner = simple(int(G, y));
% definite integral
Iinner = simple(subs(Iinner, y, d2) - subs(Iinner, y, eps));
% limit
Iinner = simple(limit(Iinner, eps, 0));

% primitive function
Iouter = simple(int(Iinner, x));
% definite integral
Iouter = simple(subs(Iouter, x, d1) - subs(Iouter, x, eps));
% limit
Iouter = simple(limit(Iouter, eps, 0));

pretty(Iouter)

%% reconstruction
R3 = sqrt(d1^2 + 2*d1*d2*c + d2^2);
a = sym('a', 'positive');
phi = sym('phi', 'real');
q = symfun(atan(a/sin(phi)+cot(phi))-atan(cot(phi)), [a, phi]);
func = (...
    d1*d2 * (3-2*log(R3)) ...
    + c*(d1^2*log(d1/R3) + d2^2*log(d2/R3))...
    - s*(d1^2*q(d2/d1,phi) + s*d2^2*q(d1/d2, phi))...
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