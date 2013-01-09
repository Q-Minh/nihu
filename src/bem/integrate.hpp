/**
 * \file integrate.hpp
 * \brief numerical integration with quadrature
 */

#ifndef INTEGRATE_HPP
#define INTEGRATE_HPP

#include <numeric>
#include "product.hpp"
#include "iterator.hpp"

/**
 * \brief metafunction to return type of quadrature weight
 * \tparam Quad quadrature class
 */
template <class Quad>
struct quad_weight_type;

/**
 * \brief metafunction to return type of quadrature location
 * \tparam Quad quadrature class
 */
template <class Quad>
struct quad_location_type;

/**
 * \brief metafunction to return type of kernel result
 * \tparam Kernel kernel class
 */
template <class Kernel>
struct kernel_result_type;

/**
 * \brief metafunction to return type of kernel argument
 * \tparam Kernel kernel class
 */
template <class Kernel>
struct kernel_arg_type;

/**
 * \brief functor performing weighted accumulation, used with std::accumulate
 * \tparam Kernel kernel class
 * \tparam Quad quadrature class
 */
template <class Kernel, class Quad>
struct WeightedAccumulator
{
	typedef typename product_type<
		typename kernel_result_type<Kernel>::type,
		typename quad_weight_type<Quad>::type
	>::type res_t;
	res_t operator () (res_t const &x, Quad const &q) 	{ return x + Kernel::eval(q.x)*q.w; }
};

template<class InputIterator, class Kernel>
struct Integral
{
	typedef typename iterator_value_type<InputIterator>::type quad_t;
	typedef typename WeightedAccumulator<Kernel, quad_t>::res_t res_t;

	static res_t eval(InputIterator begin, InputIterator end)
	{
		return std::accumulate(begin, end, res_t(), WeightedAccumulator<Kernel, quad_t>());
	}
};

#endif

