/**
 * \file math_functions.hpp
 * \brief general mathematical functions
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */
#ifndef MATH_FUNCTIONS_HPP_INCLUDED
#define MATH_FUNCTIONS_HPP_INCLUDED

/**
 * \brief \f$ sinc(x) = \sin(x) / x \f$ function
 * \tparam T type of x
 * \param [in] x the input
 * \return sinc(x)
 */
template <class T>
static T sinc(T const &x)
{
	if (std::abs(x) > 1e-3)
		return std::sin(x) / x;
	else
		return 1.0 - x*x/6.0 * (1.0 - x*x/20.0);
}


namespace bessel
{
	template <class T>
	T J0_small(T const &z)
	{
		return T();
	}

	template <class T>
	T J0_large(T const &z)
	{
		return T();
	}

	template <class T>
	T J0(T const &z)
	{
		if (abs(z) < 5.0)
			return J0_small(z);
		else
			return J0_large(z);
	}
}


#endif // MATH_FUNCTIONS_HPP_INCLUDED
