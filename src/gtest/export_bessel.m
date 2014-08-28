%EXPORT_BESSEL export Bessel function values in C arrays for testing

r = logspace(-5, 4, 15);				% modulus
phi = (0 : 1 : 360) / 180* pi;	% angle
Z = bsxfun(@times, r, exp(1i*phi'));	% complex arguments

functions = {
	'J0', @(z)besselj(0,z),
	'J1', @(z)besselj(1,z),
	'Y0', @(z)bessely(0,z),
	'Y1', @(z)bessely(1,z),
	'K0', @(z)besselk(0,z),
	'K1', @(z)besselk(1,z),
	'H01', @(z)besselh(0,1,z),
	'H11', @(z)besselh(1,1,z),
	'H02', @(z)besselh(0,2,z),
	'H12', @(z)besselh(1,2,z)
};

for i = 1 : size(functions,1)
	fun = functions{i,2};
	D = fun(Z);
    idx = ~isinf(D) & ~isnan(D) & abs(D) < 1e100 &abs(D) > 1e-100;
	name = functions {i,1};

	fid = fopen(fullfile('data', sprintf('bessel%s.dat', name)), 'wt');
	fprintf(fid, 'double bessel%s[][4] = {\n', name);
	fprintf(fid, '\t{%.10g, %.10g, %.10g, %.10g},\n', [real(Z(idx)), imag(Z(idx)), real(D(idx)), imag(D(idx))].');
	fprintf(fid, '};\n');
	fclose(fid);
end
