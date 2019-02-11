clear;

%% Define a curved line element
shapeset = ShapeSet.QuadraticLine;
coords = [
    -1, 0, 1
    -.5, 0, -.5
    ];
x = [0 .1];     % the singular point

rho0 = be_invmap(coords, x, shapeset); % and the singular point in intrinsic coordinates
zeta = rho0(2);
rho0 = rho0(1);

%% regular integral with Gaussian quadrature
[rho, w] = gaussquad(100);
[N, dN] = shapeset.eval(rho);
y = N * coords.';
jvec = dN * coords.';
jac = sqrt(dot(jvec, jvec, 2));
rvec = bsxfun(@minus, y, x);
r = sqrt(dot(rvec, rvec, 2));
G = -log(r)/2/pi;

%% expression of the singular part
[N0, dN0, ddN0] = shapeset.eval(rho0);
xs = N0 * coords.';

d0vec = xs - x;     %normal distance
d0 = norm(d0vec);
% location derivatives in the singular point
r1vec = dN0 * coords.';
r1 = sqrt(dot(r1vec, r1vec, 2));
r2vec = ddN0 * coords.';

b2 = r1^2 + dot(r2vec, d0vec, 2);
% the singular part of the integral
G0 = -log(sqrt(d0^2 + rho.^2 .* b2))/2/pi;

%%
figure;
formatfig();
plot(rho, G, rho, G0, rho, G-G0);
legend('Curved element', 'Analytically', 'difference', 'Location', 'NorthEast');
setfig('FontSize', 12, 'LineWidth', 1);
printpdf('Printed');