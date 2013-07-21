/**
 * \file product_type.hpp
 * \brief product type calculations
 * \author Peter Fiala fiala@hit.bme.hu and Peter Rucz rucz@hit.bme.hu
 */

#ifndef PRODUCT_TYPE_HPP_INCLUDED
#define PRODUCT_TYPE_HPP_INCLUDED

template<class Lhs, class Rhs>
struct product_type
{
	typedef decltype(
		(*static_cast<typename std::decay<Lhs>::type *>(nullptr))
		*
		(*static_cast<typename std::decay<Rhs>::type *>(nullptr))
	) type;
};

#endif // PRODUCT_TYPE_HPP_INCLUDED

