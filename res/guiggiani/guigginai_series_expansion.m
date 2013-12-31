function [Fm2, Fm1, beta, gamma] = guigginai_series_expansion(theta, xi0, id, X)

% location and its derivatives at the collocation point
[~, dL, ddL] = shapefun(xi0, id);
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

