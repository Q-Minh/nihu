/**
 * \file kernel_input.hpp
 * \ingroup intop
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief implementation of kernel inputs
 */

#ifndef KERNELINPUT_HPP_INCLUDED
#define KERNELINPUT_HPP_INCLUDED

#include "../tmp/bool.hpp"
#include "element.hpp"
#include "quadrature.hpp"

/** \brief traits class of kernel_inputs
 * \tparam Derived the CRTP derived class
 */
template <class Derived>
struct kernel_input_traits;

/** \brief CRTP base class of all kernel inputs */
template <class Derived>
class kernel_input_base
{
public:
	/** \brief the traits class of the Derived kernel input */
	typedef kernel_input_traits<Derived> traits_t;
	
	/** \brief the space class of the kernel input */
	typedef typename traits_t::space_t space_t;
	/** \brief the space class of the kernel input */
	typedef typename space_t::scalar_t scalar_t;
	
	/** \brief templated constructor from element and reference domain location
	 * \tparam elem_t the element type
	 */
	template <class elem_t>
	kernel_input_base(elem_t const &, typename elem_t::xi_t const &)
	{
		// check if element and kernel_input spaces are the same
		static_assert(std::is_same<space_t, typename elem_t::space_t>::value,
			"Element and kernel input spaces must match");
	}
};


#endif // KERNELINPUT_HPP_INCLUDED

