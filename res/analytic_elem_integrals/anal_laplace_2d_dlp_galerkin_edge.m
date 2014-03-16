clear all;

x = sym('x', 'positive');   % location on x elem
y = sym('y', 'positive');   % location on y elem
c = sym('c', 'real');       % cosine of angle between elements
s = sym('s', 'real');       % sine of angle between elements
d1 = sym('d1', 'positive'); % length of x elem
d2 = sym('d2', 'positive'); % length of y elem
eps = sym('eps', 'positive');   % limit variable

%% analytical derivation
rvec = [x+y*c, y*s];        % distance vector
r = sqrt(rvec * rvec.');    % scalar distance
n = [s, -c];                % normal on y elem
rdn = simple((rvec * n.') / r); % loc normal derivative

G = simple(-rdn/r/2/pi);    % Green derivative

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
a = sym('a', 'positive');
phi = sym('phi', 'real');
q = symfun(atan(a/sin(phi)+cot(phi))-atan(cot(phi)), [a, phi]);
R3 = sqrt(d1^2+d2^2+2*d1*d2*c);
func = (...
    d2*c * q(d1/d2, phi)...
    -d1 * q(d2/d1, phi)...
    - d2*s*log(R3./d2))/(2*pi);

%% difference
difference = Iouter - func;
phi = sym('phi', 'real');
simple(subs(subs(difference, s, sin(phi)), c, cos(phi)))

%% Taylor series
t0 = simple(limit(diff(subs(subs(func, s, sin(phi)), c, cos(phi)), phi, 0), phi, 0))/factorial(0);
t1 = simple(limit(diff(subs(subs(func, s, sin(phi)), c, cos(phi)), phi, 1), phi, 0))/factorial(1);
t2 = simple(limit(diff(subs(subs(func, s, sin(phi)), c, cos(phi)), phi, 2), phi, 0))/factorial(2);
t3 = simple(limit(diff(subs(subs(func, s, sin(phi)), c, cos(phi)), phi, 3), phi, 0))/factorial(3);
