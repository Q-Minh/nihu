#include "../util/math_functions.hpp"

#include <boost/math/constants/constants.hpp>

namespace NiHu
{
namespace bessel
{
	template <>
	std::complex<double> K<0>(std::complex<double> const &z)
	{
		using boost::math::double_constants::pi;
		return .5 * I * pi * (std::imag(z) < 0. ? H<0, 1>(I*z) : -H<0, 2>(-I*z));
	}

	template <>
	std::complex<double> K<1>(std::complex<double> const &z)
	{
		using boost::math::double_constants::pi;
		return -.5 * pi * (std::imag(z) < 0. ? H<1, 1>(I*z) : H<1, 2>(-I*z));
	}
}
}
