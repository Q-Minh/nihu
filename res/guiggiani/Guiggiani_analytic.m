clear;

%% symbolic variables
xi = sym('xi', 'real');
eta = sym('eta', 'real');
rho = sym('rho', 'positive');
theta = sym('theta', 'real');

%% parameters

% the geometrical shape set
L = [
    (1-xi)*(1-eta)
    (1+xi)*(1-eta)
    (1+xi)*(1+eta)
    (1-xi)*(1+eta)
    ] / 4;

% the field shape set
N = 1;

% element coordinates
X = [
    0 0 0
    1 0 0
    1.2 1 0
    0 1 0
    ];

% collocation point (x0 eta0)
xi0 = 0;
eta0 = 0;


%% first and second derivatives
% shape functions
dL = [diff(L, xi, 1), diff(L, eta, 1)];
ddL = [diff(dL(:,1), xi, 1), diff(dL(:,1), eta, 1), diff(dL(:,2), eta, 1)];
dN = [diff(N, xi, 1), diff(N, eta, 1)];

% coordinates
x = L.' * X;
dx = dL.' * X;
ddx = ddL.' * X;

% Jacobian
Jac = cross(dx(1,:), dx(2,:));

%% The Green's function
% collocation point
x0 = subs(x, {xi, eta}, {xi0, eta0});

% distance
rvec = simple(x - x0);
r = simple(sqrt(dot(rvec, rvec)));

% unit normal at collocation point
nx0 = subs(Jac, {xi, eta}, {xi0, eta0});
nx0 = nx0 / sqrt(dot(nx0, nx0));

% gradient of distance
gradr = simple(rvec / r);

% Green's function
G = N/r^3 * (dot(Jac, nx0) + 3 * dot(Jac, gradr)*dot(gradr, -nx0)) * rho;
G = subs(G, {xi, eta}, {xi0+rho*cos(theta), eta0+rho*sin(theta)});

%% Series expansion
Trig = [cos(theta), sin(theta)];

N0 = subs(N, {xi, eta}, {xi0, eta0});
N1 = subs(dN, {xi, eta}, {xi0, eta0}) * Trig.';

Avec = Trig * subs(dx, {xi, eta}, {xi0, eta0});
Bvec = [cos(theta)^2/2 cos(theta)*sin(theta), sin(theta)^2/2] * subs(ddx, {xi, eta}, {xi0, eta0});
A = sqrt(Avec * Avec.');

J0vec = subs(Jac, {xi, eta}, {xi0, eta0});
J1vec = Trig * subs([
    diff(Jac, xi, 1)
    diff(Jac, eta, 1)
    ], {xi, eta}, {xi0, eta0});

J0 = sqrt(dot(J0vec, J0vec));
Fm2 = simple(J0 * N0 / (A^3));
Fm1 = simple(((N1 * J0vec  + N0 * J1vec) / (A^3) - (3*N0*J0vec*dot(Avec, Bvec))/(A^5)) * nx0.');

%% Regularisation
Reg = simple(G - Fm2/rho^2 - Fm1/rho);

Rho = 1e-1 : 1e-1 : 1;
Theta = linspace(-pi, pi, 1e2);
[Rho, Theta] = meshgrid(Rho, Theta);

RReg = double(subs(Reg(1), {rho, theta}, {Rho, Theta}));

%% plot Taylor expansion functions and regular part
th = (-pi:1e-2:pi);
F2 = subs(Fm2, theta, th);
F1 = subs(Fm1, theta, th);

figure;
formatfig([31 8], [1 1]);

a1 = subplot(1,3,1);
plot(th, F2);
xlabel('\theta');
ylabel('F_{-2}(\theta)');

a2 = subplot(1,3,2);
plot(th, F1);
xlabel('\theta');
ylabel('F_{-1}(\theta)');

a3 = subplot(1,3,3);
surf(Theta, Rho, RReg); shading interp;
xlabel('\theta');
ylabel('\rho');
zlabel('The regular part');
