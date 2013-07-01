/**
 * \file kernel_input.hpp
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief implementation of kernel inputs
 * \details Fundamental solutions (kernels) are evaluated on two kernel inputs.
 * In the NiHu approach, the kernel defines its input that can be a single location
 * or a location with a normal. For practical purposes, the kernel input stores the Jacobian too.
 * With this approach, the kernel can make complexity estimations based on
 * the location and Jacobian of a one-point quadrature (elem center).
 */

#ifndef KERNELINPUT_HPP_INCLUDED
#define KERNELINPUT_HPP_INCLUDED

#include "element.hpp"
#include "quadrature.hpp"

/** \brief traits class of kernel_inputs */
template <class Derived>
struct kernel_input_traits;

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
	
	template <class elem_t>
	kernel_input_base(elem_t const &elem, typename elem_t::xi_t const &xi)
	{
		static_assert(std::is_same<space_t, typename elem_t::space_t>::value,
			"Element and kernel input spaces must match");
	}
};


template <class space_t>
class location;


template <class Space>
struct kernel_input_traits<location<Space> >
{
	typedef Space space_t;
};


template <class Space>
class location : public kernel_input_base<location<Space> >
{
public:
	typedef kernel_input_base<location<Space> > base_t;
	typedef typename base_t::space_t space_t;
	typedef typename space_t::location_t x_t;

	template <class elem_t>
	location(elem_t const &elem, typename elem_t::xi_t const &xi)
		: base_t(elem, xi)
	{
		m_x = elem.get_x(xi);
	}

	x_t const &get_x(void) const
	{
		return m_x;
	}

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
	x_t m_x;	/**< \brief the location */
};


template <class space_t>
class location_with_normal;


template <class Space>
struct kernel_input_traits<location_with_normal<Space> >
{
	typedef Space space_t;
};


template <class scalar_t>
class jacobian
{
public:
	jacobian()
	{
	}

	template <class elem_t>
	jacobian(elem_t const &elem, typename elem_t::xi_t const &xi)
		: m_jacobian(elem.get_normal(xi).norm())
	{
	}

	scalar_t const &get_jacobian() const
	{
		return m_jacobian;
	}

	void set_jacobian(scalar_t const &jac)
	{
		m_jacobian = jac;
	}

	template <class elem_t>
	scalar_t const &eval_on_elem(elem_t const &elem, typename elem_t::xi_t const &xi)
	{
		m_jacobian = elem.get_normal(xi).norm();
		return m_jacobian;
	}

private:
	scalar_t m_jacobian;
};


template <class Space>
class location_with_normal :
	public kernel_input_base<location_with_normal<Space> >,
	public jacobian<typename Space::scalar_t>
{
public:
	typedef kernel_input_base<location_with_normal<Space> > base_t;
	typedef typename base_t::space_t space_t;
	typedef typename space_t::location_t x_t;

	template <class elem_t>
	location_with_normal(elem_t const &elem, typename elem_t::xi_t const &xi)
		: base_t(elem, xi)
	{
		m_x = elem.get_x(xi);
		m_unit_normal = elem.get_normal(xi);
		this->set_jacobian(m_unit_normal.norm());
		m_unit_normal /= this->get_jacobian();
	}

	x_t const &get_x(void) const
	{
		return m_x;
	}

	x_t const &get_unit_normal(void) const
	{
		return m_unit_normal;
	}

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
	x_t m_x;	/**< \brief the location */
	x_t m_unit_normal;	/**< \brief the unit normal vector */
};



template <class KernelInput>
class weighted_kernel_input_impl :
	public KernelInput,
	public jacobian<typename KernelInput::scalar_t>
{
public:
	typedef typename KernelInput::scalar_t scalar_t;

	template<class elem_t>
	weighted_kernel_input_impl(elem_t const &elem, typename elem_t::xi_t const &xi) :
		KernelInput(elem, xi),
		jacobian<scalar_t>(elem, xi)
	{
	}
};


template <bool HasJacobian, class KernelInput>
struct weighted_kernel_input_selector
{
	typedef KernelInput type;
};


template <class KernelInput>
struct weighted_kernel_input_selector<false, KernelInput>
{
	typedef weighted_kernel_input_impl<KernelInput> type;
};


template <class KernelInput>
struct weighted_kernel_input : weighted_kernel_input_selector<
	std::is_base_of<
		jacobian<typename KernelInput::scalar_t>,
		KernelInput
	>::value,
	KernelInput
> {};


#endif // KERNELINPUT_HPP_INCLUDED

