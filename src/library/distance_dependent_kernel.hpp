/** \file distance_dependent_kernel.hpp
 * \brief implementation of class distance_dependent_kernel
 * \author Peter Fiala
 */

#ifndef DISTANCE_DEPENDENT_KERNEL_HPP_INCLUDED
#define DISTANCE_DEPENDENT_KERNEL_HPP_INCLUDED

#include "../util/crtp_base.hpp"

namespace NiHu
{

template <class Derived>
class distance_dependent_kernel;

namespace distance_dependent_kernel_traits_ns
{
	template <class Derived>
	struct space;

	template <class Derived>
	struct result;

	template <class Derived>
	struct quadrature_family;

	template <class Derived>
	struct is_singular;

	template <class DerivedSpace>
	struct singular_core;

	template <class Derived>
	struct singular_quadrature_order;

	template <class Derived>
	struct far_field_behaviour;

	template <class Derived>
	struct singularity_type;
} // end of namespace distance_dependent_kernel_traits_ns

template <class Derived>
class distance_dependent_kernel
{
public:
	NIHU_CRTP_HELPERS
};

} // end of namespace NiHu

#endif // DISTANCE_DEPENDENT_KERNEL_HPP_INCLUDED
