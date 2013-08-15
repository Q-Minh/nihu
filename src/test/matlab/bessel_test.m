clear;
z = linspace(-10, 10, 1001)';
J0 = bessel(complex(z));

plot(z, [real(J0), imag(J0), real(besselh(0,z)), imag(besselh(0,z))]);