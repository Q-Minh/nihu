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

template <class T>
static T besselj0_small(T const &z)
{
	T res = 1.0, q(z*z/4.0);
	for (unsigned n = 10; n > 0; --n)
		res = 1.0 - q/(n*n)*res;
	return res;
}

template <class T>
static T besselj1_small(T const &z)
{
	T res = 1.0, q(z*z/4.0);
	for (unsigned n = 10; n > 0; --n)
		res = 1.0 - q/(n*(n+1))*res;
	return res * z/2.0;
}

#endif // MATH_FUNCTIONS_HPP_INCLUDED
