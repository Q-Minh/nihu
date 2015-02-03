function g = elastodynamic_pulsating_sphere(R0, R, mu, rho, nu, freq)

om = 2*pi*freq;
cS = sqrt(mu/rho);
a = sqrt((1-2*nu)./(2*(1-nu)));
kS = om/cS;
kP = a * kS;

g = R0^3./R.^2 / (4*mu) .* (1+1i*kP*R) / (1 + 1i*kP*R0 - (kS*R0/2)^2) .* exp(-1i*(R-R0)*kP);

end
