#include <library/log_gauss_quadrature.hpp>
#include <iostream>

double func(double x, unsigned n)
{
	double res = -std::log(x);
	for (unsigned k = 0; k < n; ++k)
		res *= x;
	return res;
}

int main(void)
{
	unsigned N = 5;
	auto quadrature = log_gauss_impl<double>(N);

	for (unsigned n = 0; n < 2 * N + 5; ++n)
	{
		double res = 0.0;
		for (unsigned i = 0; i < quadrature.rows(); ++i)
			res += func(quadrature(i, 0), n) * quadrature(i, 1);
		double anal = 1.0 / (n + 1) / (n + 1);
		std::cout << n << '\t' << std::log10(std::abs(res / anal - 1.0)) << std::endl;
	}

	return 0;
}
