function [I, Isurf, Ilin1, Ilin2] = Guiggiani(n)

X = [
    -1 -1 0
    1 -1 0
    1 1 0
    -1 1 0
    ];
xi0 = [0 0];

%% surface integral
% create surface Duffy quadrature
[xi, w] = DuffyQuad(xi0, n);
% the same in polar
[theta, rho] = cart2pol(xi(:,1)-xi0(1), xi(:,2)-xi0(2));
w = w./ rho;

[Fm2, Fm1] = series_expansion(theta, xi0, X);
Series = bsxfun(@times, Fm2, 1./rho.^2) + bsxfun(@times, Fm1, 1./rho);

% the hypersingular function

% location and normal at the collocational point
[L, dL] = shapefun(xi0, 24);
x0 = L * X;
dx0_xi = dL(:,:,1) * X;
dx0_eta = dL(:,:,2) * X;
nx0 = cross(dx0_xi, dx0_eta);
nx0 = nx0 ./ norm(nx0);

[L, dL] = shapefun(xi, 24);
% location and its derivatives over the element
x = L * X;
dx_xi = dL(:,:,1) * X;
dx_eta = dL(:,:,2) * X;
% Jacobian vector
J = cross(dx_xi, dx_eta);
% distance, scalar distance and grad r
rvec = bsxfun(@minus, x, x0);
r = sqrt(dot(rvec, rvec, 2));
gradr = bsxfun(@times, rvec, 1./r);
% shape function values
N = shapefun(xi, 21);
% G'' * N * rho
F = ((J*nx0.')-3*(dot(gradr,J,2)).*(-gradr*nx0.')).*rho./ (4*pi*r.^3);
F = bsxfun(@times, F, N);

Isurf = w' * (F - Series);

%% line integrals
Ilin1 = 0;
Ilin2 = 0;
coords = [
    -1 -1
    1 -1
    1 1
    -1 1
    ];
theta_lim = atan2(coords(:,2)-xi0(2), coords(:,1)-xi0(1));
for l = 1 : 4
    t1 = theta_lim(l);
    t2 = theta_lim(mod(l+1-1, 4)+1);
    if (t2 < t1)
        t2 = t2 + 2*pi;
    end
    [theta, w] = gaussquad(n, t1, t2);

    switch l
        case 1
            rho_lim = (-1 - xi0(2)) ./ sin(theta);
        case 2
            rho_lim = (+1 - xi0(1)) ./ cos(theta);
        case 3
            rho_lim = (+1 - xi0(2)) ./ sin(theta);
        case 4
            rho_lim = (-1 - xi0(1)) ./ cos(theta);
    end
    
    [Fm2, Fm1, beta, gamma] = series_expansion(theta, xi0, X);

    Ilin1 = Ilin1 + w' * (bsxfun(@times, Fm1, log(abs(rho_lim ./ beta))));
    Ilin2 = Ilin2 + w' * (bsxfun(@times, Fm2, gamma./beta.^2 + 1./rho_lim));
end

%% The total contribution
I = Isurf + Ilin1 - Ilin2;

end

function [Fm2, Fm1, beta, gamma] = series_expansion(theta, xi0, X)

% location and its derivatives at the collocation point
[~, dL, ddL] = shapefun(xi0, 24);
dx0_xi = dL(:,:,1) * X;
dx0_eta = dL(:,:,2) * X;
ddx0_xixi = ddL(:,:,1) * X;
ddx0_xieta = ddL(:,:,2) * X;
ddx0_etaeta = ddL(:,:,3) * X;

% jacobian and its derivatives at the collocation point
J0vec = cross(dx0_xi, dx0_eta);
J0 = norm(J0vec);
nx0 = J0vec / J0;
J1vec = cos(theta) * (cross(ddx0_xixi, dx0_eta) + cross(dx0_xi, ddx0_xieta)) +...
    sin(theta) * (cross(ddx0_xieta, dx0_eta) + cross(dx0_eta, ddx0_etaeta));

Avec = cos(theta) * dx0_xi + sin(theta) * dx0_eta;
Bvec = cos(theta).^2/2 * ddx0_xixi +...
    sin(theta).*cos(theta) * ddx0_xieta +...
    sin(theta).^2/2 * ddx0_etaeta;
A = sqrt(dot(Avec, Avec, 2));

% field shape functions
[N0, dN] = shapefun(xi0, 21);
N1 = cos(theta) * dN(:,:,1) + sin(theta) * dN(:,:,2);

Fm2 = (1./(4*pi*A.^3)) * (J0 .* N0);
Fm1 = bsxfun(@times, 1./(4*pi*A.^3), (J0vec*nx0.') .* N1) +...
    bsxfun(@times, 1./(4*pi*A.^3), (J1vec*nx0.') * N0) -...
    3*bsxfun(@times, 1./(4*pi*A.^5), dot(Avec, Bvec,2) * (J0vec*nx0.') * N0);

if (nargout > 2)
    beta = 1./A;
    gamma = -dot(Avec, Bvec, 2)./A.^4;
end
end

function [Xi, W] = DuffyQuad(xi0, n)

order = 2*n-1;
[xi, w] = gaussquad2(order);
c = [
    -1 -1
    1 -1
    1 1
    -1 1
    ];

[N, dN] = shapefun(xi, 24);

Xi = [];
W = [];

for i = 1 : 4
    C(1,:) = xi0;
    C(2,:) = xi0;
    C(3,:) = c(i,:);
    C(4,:) = c(mod(i+1-1, 4)+1,:);
    
    Xi = [
        Xi
        N * C
        ];
    dx1 = dN(:,:,1) * C;
    dx2 = dN(:,:,2) * C;
    ww = dx1(:,1).*dx2(:,2) - dx2(:,1).*dx1(:,2);
    W = [
        W
        w.*ww
        ];
end

end
