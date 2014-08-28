clear;

r = logspace(-4, 2, 7);
phi = (-5 : 15 : 360) / 180* pi;
[R, PHI] = meshgrid(r, phi);
Z = R .* exp(1i*PHI);

D = bessely(1, Z);

fid = fopen('data/besselY1.dat', 'wt');
fprintf(fid, 'double besselY1[][4] = {\n');
fprintf(fid, '\t{%.10g, %.10g, %.10g, %.10g},\n', [real(Z(:)), imag(Z(:)), real(D(:)), imag(D(:))].');
fprintf(fid, '};\n');
fclose(fid);
