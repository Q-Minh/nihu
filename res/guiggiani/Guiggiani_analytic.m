clear;

%% symbolic variables
xi = sym('xi', 'real');
eta = sym('eta', 'real');
rho = sym('rho', 'positive');
theta = sym('theta', 'real');
Lx = sym('Lx', 'positive');
% Ly = sym('Ly', 'positive');
Ly = Lx;

%% parameters

% the geometrical shape set
L = [
    (1-xi)*(1-eta)
    (1+xi)*(1-eta)
    (1+xi)*(1+eta)
    (1-xi)*(1+eta)
    ] / 4;

% the field shape set
N = [
    (1.0 - sqrt(3.0)*xi) * (1.0 - sqrt(3.0)*eta)
    (1.0 + sqrt(3.0)*xi) * (1.0 - sqrt(3.0)*eta)
    (1.0 - sqrt(3.0)*xi) * (1.0 + sqrt(3.0)*eta)
    (1.0 + sqrt(3.0)*xi) * (1.0 + sqrt(3.0)*eta)
    ] / 4;

% element coordinates
X = [
    0 0 0
    Lx 0 0
    Lx Ly 0
    0 Ly 0
    ];

% collocation point (xi0 eta0)
xi0 = -1/sqrt(3);
eta0 = -1/sqrt(3);

%% first and second derivatives
% shape functions
dL = [diff(L, xi, 1), diff(L, eta, 1)];
ddL = [diff(dL(:,1), xi, 1), diff(dL(:,1), eta, 1), diff(dL(:,2), eta, 1)];
dN = [diff(N, xi, 1), diff(N, eta, 1)];

% coordinates
x = L.' * X;
dx = dL.' * X;
ddx = ddL.' * X;

% Jacobian vector
Jvec = cross(dx(1,:), dx(2,:));

%% The Green's function
% collocation point
x0 = subs(x, {xi, eta}, {xi0, eta0});

% distance
rvec = simplify(x - x0);
r = simplify(sqrt(dot(rvec, rvec)));

% unit normal at collocation point
J0vec = subs(Jvec, {xi, eta}, {xi0, eta0});
nx0 = J0vec / sqrt(dot(J0vec, J0vec));

% gradient of distance
gradr = simplify(rvec / r);

% integrand
F = N/r^3 * (dot(Jvec, nx0) + 3 * dot(Jvec, gradr)*dot(gradr, -nx0));
% integrand
F = subs(F, {xi, eta}, {xi0+rho*cos(theta), eta0+rho*sin(theta)}) * rho;

%% Series expansion
Trig = [cos(theta), sin(theta)];

N0 = subs(N, {xi, eta}, {xi0, eta0});
N1 = subs(dN, {xi, eta}, {xi0, eta0}) * Trig.';

Avec = Trig * subs(dx, {xi, eta}, {xi0, eta0});
Bvec = [cos(theta)^2/2 cos(theta)*sin(theta), sin(theta)^2/2] * subs(ddx, {xi, eta}, {xi0, eta0});
A = sqrt(Avec * Avec.');

J1vec = Trig * subs([
    diff(Jvec, xi, 1)
    diff(Jvec, eta, 1)
    ], {xi, eta}, {xi0, eta0});

J0 = sqrt(dot(J0vec, J0vec));
Fm2 = simplify(J0 * N0 / (A^3));
Fm1 = simplify(((N1 * J0vec  + N0 * J1vec) / (A^3) - (3*N0*J0vec*dot(Avec, Bvec))/(A^5)) * nx0.');

%% Regularisation
FReg = simplify(F - Fm2/rho^2 - Fm1/rho);


%% numerical integral
Im2 = 0;
Im1 = 0;
Ir = 0;
[theta_lim, theta0, rho0] = plane_elem_helper_mid([-1 -1; 1 -1; 1 1; -1 1], [xi0 eta0]);
order = 100;
for i = 1 : 4
    t1 = theta_lim(i);
    t2 = theta_lim(mod(i,4)+1);
    if (abs(t2-t1) > pi)
        if (t2 < t1)
            t2 = t2 + 2*pi;
        else
            t1 = t1 + 2*pi;
        end
    end
    
    [th, w_th] = gaussquad(order, t1, t2);
    th = th(:).';
    rholim = rho0(i)./cos(th-theta0(i));
    
    im2 = double(subs(-Fm2/rho, {rho, theta, Lx}, {rholim, th, 1}));
    im1 = double(subs(Fm1*log(abs(rho)), {rho, theta, Lx}, {rholim, th, 1}));
    ir = double(subs(int(FReg, rho), {rho, theta, Lx}, {rholim, th, 1}));
    
    Im2 = Im2 + im2 * w_th;
    Im1 = Im1 + im1 * w_th;
    Ir = Ir + ir * w_th;
end

I = Im2 + Im1 + Ir;

