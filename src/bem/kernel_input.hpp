/**
 * \file kernel_input.hpp
 * \ingroup kernel
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief implementation of kernel inputs
 */

#ifndef KERNELINPUT_HPP_INCLUDED
#define KERNELINPUT_HPP_INCLUDED

#include "element.hpp"
#include "quadrature.hpp"

/** \brief traits class of kernel_inputs
 *
 * \snippet custom_kernel.cpp Kernel input traits
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
	 * \param [in] elem the element
	 * \param [in] xi point in the reference domain
	 */
	template <class elem_t>
	kernel_input_base(elem_t const &elem, typename elem_t::xi_t const &xi)
	{
		// check if element and kernel_input spaces are the same
		static_assert(std::is_same<space_t, typename elem_t::space_t>::value,
			"Element and kernel input spaces must match");
	}
};


/** \brief class to store the jacobian of the element's coordinate transform
 * \tparam scalar_t the scalar type of the Jacobian
 */
template <class scalar_t>
class jacobian
{
public:
	/** \brief default constructor */
	jacobian()
	{
	}

	/** \brief constructor from element and reference domain coordinate
	 * \tparam elem_t the element type
	 * \param [in] elem reference to an element
	 * \param [in] xi reference to a location in the xi domain
	 */
	template <class elem_t>
	jacobian(elem_t const &elem, typename elem_t::xi_t const &xi)
		: m_jacobian(elem.get_normal(xi).norm())
	{
	}

	/** \brief return Jacobian
	 * \return the stored Jacobian
	 */
	scalar_t const &get_jacobian() const
	{
		return m_jacobian;
	}

	/** \brief set Jacobian
	 * \param [in] jac the new Jacobian
	 */
	void set_jacobian(scalar_t const &jac)
	{
		m_jacobian = jac;
	}

	/** \brief evaluate Jacobian on an element
	 * \tparam elem_t the element type
	 * \param [in] elem reference to an element
	 * \param [in] xi reference to a location in the xi domain
	 */
	template <class elem_t>
	scalar_t const &eval_on_elem(elem_t const &elem, typename elem_t::xi_t const &xi)
	{
		m_jacobian = elem.get_normal(xi).norm();
		return m_jacobian;
	}

private:
	/** \brief the stored Jacobian value */
	scalar_t m_jacobian;
};


/** \brief metafunction implementation to assign a weighted kernel input to a kernel input */
template <class KernelInput>
class weighted_kernel_input_impl :
	public KernelInput,
	public jacobian<typename KernelInput::scalar_t>
{
public:
	/** \brief constructor from element and a reference domain variable
	 * \tparam elem_t the element type
	 * \param [in] elem the element reference
	 * \param [in] xi the reference domain variable
	 */
	template<class elem_t>
	weighted_kernel_input_impl(elem_t const &elem, typename elem_t::xi_t const &xi) :
		KernelInput(elem, xi),
		jacobian<typename KernelInput::scalar_t>(elem, xi)
	{
	}
};

/** \brief metafunction to assign a weighted kernel input to a kernel input */
template <class KernelInput>
struct weighted_kernel_input : tmp::if_<
	typename std::is_base_of<
		jacobian<typename KernelInput::scalar_t>,
		KernelInput
	>::type,
	KernelInput,
	weighted_kernel_input_impl<KernelInput>
> {};

#endif // KERNELINPUT_HPP_INCLUDED
