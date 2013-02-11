#ifndef DESCRIPTOR_HPP_INCLUDED
#define DESCRIPTOR_HPP_INCLUDED

#include <type_traits>

#include "element.hpp"
#include "quadrature.hpp"

/**
 * \brief location class containing a locatin and a jacobian
 * \tparam xType type of the location
 */
template <class xType>
class location
{
public:
	typedef xType x_t;	/**< \brief template argument as nested type */
	typedef typename x_t::Scalar scalar_t; /**< \brief the scalar type of the descriptor */

	/**
	 * \brief default constructor needed for elem_accelerator
	 */
	location()
	{
	}

	/**
	 * \brief constructor from an element and a quadrature point
	 * \tparam elem_t
	 * \param elem the element
	 * \param q a quadrature elem
	 */
	template <class elem_t>
	location(elem_t const &elem, quadrature_elem<typename elem_t::domain_t> const &q)
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
	 * \brief return the location
	 * \return the location
	 */
	x_t const &get_x(void) const
	{
		return m_x;
	}

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
	x_t m_x;	/**< \brief the stored location */
	scalar_t m_jacobian; /**< \brief the stored jacobian */
};


/**
 * \brief location class containing a location, a normal and a jacobian
 * \tparam xType type of the location
 */
template <class xType>
class location_with_normal : public location<xType>
{
public:
	typedef xType x_t; 	/**< \brief template argument as nested type */
	typedef location<xType> base; 	/**< \brief the base class */

	/**
	 * \brief default constructor needed for elem_accelerator
	 */
	location_with_normal()
	{
	}

	/**
	 * \brief constructor from an element and a quadrature point
	 * \tparam elem_t
	 * \param elem the element
	 * \param q a quadrature elem
	 */
	template <class elem_t>
	location_with_normal(elem_t const &elem, quadrature_elem<typename elem_t::domain_t> const &q)
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


#endif

