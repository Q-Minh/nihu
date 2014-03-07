clear all;

y = sym('y', 'real');
x = sym('x', 'real');
d = sym('d', 'positive');
eps = sym('eps', 'positive');

%% green's function
r = abs(y-x);
G = -log(r)/2/pi;

%% analytical derivation for constant shape functions
I0 = int(limit(int(G, y, 0, x-eps) + int(G, y, x+eps, d), eps, 0), x, 0, d);
pretty(I0)

%% analytical derivation for linear shape functions
Nx = [1-x/d, x/d];
Ny = [1-y/d, y/d];
I = int(Nx.' * limit(int(G * Ny, y, 0, x-eps) + int(G * Ny, y, x+eps, d), eps, 0), x, 0, d);
pretty(I)

%% green's function
r = abs(y-x);
G = 1/r^2/2/pi;

%% analytical derivation for constant shape functions
I0 = int(int(G, y, 0, x-eps) + int(G, y, x+eps, d), x, 0, d);
pretty(I0)

%% analytical derivation for linear shape functions
Nx = [1-x/d, x/d];
Ny = [1-y/d, y/d];
I = int(Nx.' * limit(int(G * Ny, y, 0, x-eps) + int(G * Ny, y, x+eps, d), eps, 0), x, 0, d);
pretty(I)

