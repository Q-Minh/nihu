function HyperLaplaceTest

X = [
    0 0 0
    1 0 0
    1 1 0
    0 1.2 .2
    ];
xi0 = [.5 .5];
n = 7;

% create surface Duffy quadrature
[xi, w] = DuffyQuad(xi0, n);
% the same in polar
[theta, rho] = cart2pol(xi(:,1)-xi0(1), xi(:,2)-xi0(2));

% location and its derivatives at the collocation point
[L, dL, ddL] = shapefun(xi0, 24);
x0 = L * X;
dx0_xi = dL(:,:,1) * X;
dx0_eta = dL(:,:,2) * X;
ddx0_xixi = ddL(:,:,1) * X;
ddx0_xieta = ddL(:,:,2) * X;
ddx0_etaeta = ddL(:,:,3) * X;

Avec = cos(theta) * dx0_xi + sin(theta) * dx0_eta;
Bvec = cos(theta).^2/2 * ddx0_xixi +...
    sin(theta).*cos(theta) * ddx0_xieta +...
    sin(theta).^2/2 * ddx0_etaeta;
A = sqrt(dot(Avec, Avec, 2));

% jacobian and its derivatives at the collocation point
J0vec = cross(dx0_xi, dx0_eta);
J0 = norm(J0vec);
nx0 = J0vec / J0;
J1vec = cos(theta) * (cross(ddx0_xixi, dx0_eta) + cross(dx0_xi, ddx0_xieta)) +...
    sin(theta) * (cross(ddx0_xieta, dx0_eta) + cross(dx0_eta, ddx0_etaeta));

% field shape functions
[N0, dN] = shapefun(xi0, 24);
N1 = cos(theta) * dN(:,:,1) + sin(theta) * dN(:,:,2);

Fm2 = (1./(4*pi*A.^3)) * (J0 .* N0);
Fm1 = bsxfun(@times, 1./(4*pi*A.^3), (J0vec*nx0.') .* N1) +...
    bsxfun(@times, 1./(4*pi*A.^3), (J1vec*nx0.') * N0) -...
    3*bsxfun(@times, 1./(4*pi*A.^5), dot(Avec, Bvec,2) * (J0vec*nx0.') * N0);
Series = bsxfun(@times, Fm2, 1./rho.^2) + bsxfun(@times, Fm1, 1./rho);


% the hypersingular function
[L, dL] = shapefun(xi, 24);
x = L * X;
dx_xi = dL(:,:,1) * X;
dx_eta = dL(:,:,2) * X;
J = cross(dx_xi, dx_eta);
n = bsxfun(@times, J, 1./sqrt(dot(J,J,2)));
rvec = bsxfun(@minus, x, x0);
r = sqrt(dot(rvec, rvec, 2));
gradr = bsxfun(@times, rvec, 1./r);
N = shapefun(xi, 24);
F = ((J*nx0.')-3*(dot(gradr, J,2)).*(-gradr*nx0.')).*rho./ (4*pi*r.^3);
F = bsxfun(@times, F, N);

Regular = F - Series;

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
