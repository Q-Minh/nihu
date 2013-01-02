/**
 * \file integrate.hpp
 * \brief numerical integration with quadrature
 */

#ifndef INTEGRATE_HPP
#define INTEGRATE_HPP

#include <numeric>

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
 * \brief metafunction to return the type of a product
 * \tparam A left hand side term type
 * \tparam B right hand side term type
 */
template <class A, class B>
struct product_type;

/**
 * \brief functor performing weighted accumulation, used with std::accumulate
 * \tparam Kernel kernel class
 * \tparam Quad quadrature class
 */
template <class Kernel, class Quad>
struct IntegFunctor
{
	typedef typename product_type<
		typename kernel_result_type<Kernel>::type,
		typename quad_weight_type<Quad>::type
	>::type res_t;
	res_t operator () (res_t const &x, Quad const &q) 	{ return x + Kernel::apply(q.x)*q.w; }
};

#endif

