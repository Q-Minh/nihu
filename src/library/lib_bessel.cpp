#include "../util/math_functions.hpp"


namespace NiHu
{
namespace bessel
{
	template <>
	std::complex<double> K<0>(std::complex<double> const &z)
	{
		return .5*I*M_PI* (std::imag(z) < 0. ? H<0, 1>(I*z) : -H<0, 2>(-I*z));
	}

	template <>
	std::complex<double> K<1>(std::complex<double> const &z)
	{
		return -.5*M_PI * (std::imag(z) < 0. ? H<1, 1>(I*z) : H<1, 2>(-I*z));
	}
}
}
