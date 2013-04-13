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

/**
 * \brief class representing a kernel input consisting of a location and a jacobian
 * \details A location is constructed from and element and a quadrature point.
 * The jacobian incorporates the quadrature weight.
 * \tparam xType location type of the element
 */
template <class xType>
class location
{
protected:
	/**
	 * \brief default constructor only needed to avoid constructor calls in derived classes
	 */
	location()
	{
	}

public:
	typedef xType x_t;	/**< \brief template parameter as nested type */
	typedef typename x_t::Scalar scalar_t;	/**< \brief the scalar type of the location */

	/**
	 * \brief constructor from element and quadrature point
	 * \tparam elem_t the element type
	 * \param elem an element
	 * \param q a quadrature point
	 */
	template <class elem_t>
	location(elem_t const &elem, quadrature_elem<typename elem_t::domain_t::xi_t, typename elem_t::domain_t::scalar_t> const &q)
	{
		static_assert(std::is_same<x_t, typename elem_t::x_t>::value,
			"Element and descriptor location types must match");
		typename elem_t::xi_t xi = q.get_xi();
		m_x = elem.get_x(xi);
		m_jacobian = elem.get_normal(xi).norm() * q.get_w();
	}

	/**
	 * \brief return the jacobian
	 * \return the jacobian
	 */
	scalar_t const &get_jacobian(void) const
	{
		return m_jacobian;
	}

	/**
	 * \brief return location
	 * \return location
	 */
	x_t const &get_x(void) const
	{
		return m_x;
	}

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
	x_t m_x;	/**< \brief the location */
	scalar_t m_jacobian;	/**< \brief the jacobian */
};


/**
 * \brief class representing a kernel input consisting of a location, normal vector and a jacobian
 * \tparam xType location type of the element
 * \details A location is constructed from and element and a quadrature point.
 * The jacobian incorporates the quadrature weight.
 */
template <class xType>
class location_with_normal : public location<xType>
{
public:
	typedef xType x_t;	 /**< \brief template  parameter as nested type */
	typedef location<xType> base;	/**< \brief the base class */

	/**
	 * \brief constructor from element and quadrature point
	 * \tparam elem_t the element type
	 * \param elem an element
	 * \param q a quadrature point
	 */
	template <class elem_t>
	location_with_normal(elem_t const &elem, quadrature_elem<typename elem_t::domain_t::xi_t, typename elem_t::domain_t::scalar_t> const &q)
	{
		static_assert(std::is_same<x_t, typename elem_t::x_t>::value,
			"Element and descriptor location types must match");
		typename elem_t::xi_t xi = q.get_xi();
		base::m_x = elem.get_x(xi);
		m_normal = elem.get_normal(q.get_xi());
		base::m_jacobian = m_normal.norm();
		m_normal /= base::m_jacobian;
		base::m_jacobian *= q.get_w();
	}

	/**
	 * \brief return the normal
	 * \return the unit normal vector
	 */
	x_t const &get_normal(void) const
	{
		return m_normal;
	}

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
	x_t m_normal;	/**< \brief the stored unit normal vector */
};


#endif // KERNELINPUT_HPP_INCLUDED

