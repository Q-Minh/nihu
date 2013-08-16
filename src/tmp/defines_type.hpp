/**
 * \file defines_type.hpp
 * \brief metafunction determining typedefs
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef DEFINES_TYPE_HPP_INCLUDED
#define DEFINES_TYPE_HPP_INCLUDED

#include <type_traits>

namespace internal
{
	/** \brief metafunction converting anything to void */
	template <class C>
	struct voidize { typedef void type; };
}

#define DEFINES_TYPE(TYPE)								\
template <class C, class = void>						\
struct defines_##TYPE : std::false_type {};				\
template <class C>										\
														\
struct defines_##TYPE <									\
	C,													\
	typename internal::voidize<typename C::TYPE>::type	\
> : std::true_type { };

#endif // DEFINES_TYPE_HPP_INCLUDED
