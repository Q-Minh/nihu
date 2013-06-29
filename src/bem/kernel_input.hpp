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
	
	/** \brief metafunction assigning a quadrature element to an element */
	template <class element>
	struct quadr_elem
	{
		typedef quadrature_elem<
			typename element::domain_t::xi_t,
			typename element::domain_t::scalar_t
		> type;
	};
	
	template <class elem_t>
	kernel_input_base(elem_t const &elem, typename quadr_elem<elem_t>::type const &qe)
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
	typedef typename space_t::scalar_t scalar_t;

	template <class elem_t, class quadrature_elem_t>
	location(elem_t const &elem, quadrature_elem_t const &q)
		: base_t(elem, q)
	{
		typename elem_t::xi_t xi = q.get_xi();
		m_x = elem.get_x(xi);
		m_jacobian = elem.get_normal(xi).norm() * q.get_w();
	}

	scalar_t const &get_jacobian(void) const
	{
		return m_jacobian;
	}

	x_t const &get_x(void) const
	{
		return m_x;
	}

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
	x_t m_x;	/**< \brief the location */
	scalar_t m_jacobian;	/**< \brief the jacobian */
};


template <class space_t>
class location_with_normal;


template <class Space>
struct kernel_input_traits<location_with_normal<Space> >
{
	typedef Space space_t;
};



template <class Space>
class location_with_normal : public kernel_input_base<location_with_normal<Space> >
{
public:
	typedef kernel_input_base<location_with_normal<Space> > base_t;
	typedef typename base_t::space_t space_t;
	typedef typename space_t::location_t x_t;
	typedef typename space_t::scalar_t scalar_t;

	template <class elem_t, class quadrature_elem_t>
	location_with_normal(elem_t const &elem, quadrature_elem_t const &q)
		: base_t(elem, q)
	{
		typename elem_t::xi_t xi = q.get_xi();
		m_x = elem.get_x(xi);
		m_normal = elem.get_normal(q.get_xi());
		m_jacobian = m_normal.norm();
		m_normal /= m_jacobian;
		m_jacobian *= q.get_w();
	}

	scalar_t const &get_jacobian(void) const
	{
		return m_jacobian;
	}

	x_t const &get_x(void) const
	{
		return m_x;
	}

	x_t const &get_normal(void) const
	{
		return m_normal;
	}

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
	x_t m_x;	/**< \brief the location */
	scalar_t m_jacobian;	/**< \brief the jacobian */
	x_t m_normal;	/**< \brief the stored unit normal vector */
};


#endif // KERNELINPUT_HPP_INCLUDED

