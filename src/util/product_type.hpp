/**
 * \file product_type.hpp
 * \brief product type calculations
 * \author Peter Fiala fiala@hit.bme.hu and Peter Rucz rucz@hit.bme.hu
 */

#ifndef PRODUCT_TYPE_HPP_INCLUDED
#define PRODUCT_TYPE_HPP_INCLUDED

#include "plain_type.hpp"

/** \brief metafunction returning the product type of two classes
 * \tparam Lhs the left hand side expression type
 * \tparam Rhs the right hand side expression type
 */
template<class Lhs, class Rhs, bool isBothEigen = is_eigen<Lhs>::value && is_eigen<Rhs>::value>
struct product_type
{
	/** \brief the return type computed by decltype */
	typedef decltype(
		(*static_cast<typename std::decay<Lhs>::type *>(nullptr))
		*
		(*static_cast<typename std::decay<Rhs>::type *>(nullptr))
	) type;
};


template<class Lhs, class Rhs>
struct product_type<Lhs, Rhs, true>
{
	typedef typename Eigen::ProductReturnType<Lhs, Rhs>::Type type;
};

#endif // PRODUCT_TYPE_HPP_INCLUDED

