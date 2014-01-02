function [I, Isurf, Ilin1, Ilin2] = Guiggiani(n, id, X, xi0)
%GUIGGIANI  Guiggiani's method to compute hypersingular integrals
%   I = GUIGGIANI(n, X, xi0) computes the hypersingular itnegral with
%   Guiggiani's method, where the element is described by the corner
%   coordinates X, the collocation point is located at xi0, and a n-point
%   Gaussian quadrature is used.

%% surface integral
% create surface Duffy quadrature
[xi, w] = DuffyQuad(n, id, xi0);
% [xi, w] = gaussquad2(2*n-1, 4);
% the same in polar
[theta, rho] = cart2pol(xi(:,1)-xi0(1), xi(:,2)-xi0(2));
% decompensate null Jacobian
w = w./ rho;

[Fm2, Fm1] = guiggiani_series_expansion(theta, xi0, id, X);
Series = bsxfun(@times, Fm2, 1./rho.^2) + bsxfun(@times, Fm1, 1./rho);

% the hypersingular function

% location and normal at the collocational point
[L, dL] = shapefun(xi0, id);
x0 = L * X;
dx0_xi = dL(:,:,1) * X;
dx0_eta = dL(:,:,2) * X;
nx0 = cross(dx0_xi, dx0_eta);
nx0 = nx0 ./ norm(nx0);

[L, dL] = shapefun(xi, id);
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
F = ((J*nx0.')+3*(dot(gradr,J,2)).*(-gradr*nx0.')).*rho./ (4*pi*r.^3);
F = bsxfun(@times, F, N);

Isurf = w' * (F - Series);

%% line integrals
Ilin1 = 0;
Ilin2 = 0;
coords = reference_domain_corners(id);
theta_lim = atan2(coords(:,2)-xi0(2), coords(:,1)-xi0(1));
nC = mod(id,10);
for l = 1 : nC
    t1 = theta_lim(l);
    t2 = theta_lim(mod(l+1-1, nC)+1);
    if (t2 < t1)
        t2 = t2 + 2*pi;
    end
    [theta, w] = gaussquad(n, t1, t2);
    [Fm2, Fm1, beta, gamma] = guiggiani_series_expansion(theta, xi0, id, X);

    rho_lim = distance_to_reference_boundary(xi0, theta, id);
    Ilin1 = Ilin1 + w' * (bsxfun(@times, Fm1, log(abs(rho_lim ./ beta))));
    Ilin2 = Ilin2 + w' * (bsxfun(@times, Fm2, gamma./beta.^2 + 1./rho_lim));
end

%% The total contribution
I = Isurf + Ilin1 - Ilin2;

end

