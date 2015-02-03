function g = elastodynamic_green(rvec, mu, rho, nu, freq)

om = 2*pi*freq;
cS = sqrt(mu/rho);
a = sqrt((1-2*nu)./(2*(1-nu)));
r = norm(rvec);
gamma = rvec/r;
kSr = om/cS*r;
kPr = a * kSr;

psi = exp(-1i*kPr) * a^2 .* (1i/kPr + 1./kPr.^2) ...
    + exp(-1i*kSr) .* (1 - 1i./kSr - 1./kSr.^2);
chi = exp(-1i*kPr) * a^2 .* (1-3*1i./kPr - 3./kPr.^2) ...
    - exp(-1i*kSr) .* (1 - 3*1i./kSr - 3./kSr.^2);

g = (psi * eye(3) + chi * (gamma(:) * gamma(:)')) / (4*pi*r*mu);

end
