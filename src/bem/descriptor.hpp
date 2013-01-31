#ifndef DESCRIPTOR_HPP_INCLUDED
#define DESCRIPTOR_HPP_INCLUDED

#include <type_traits>

#include "element.hpp"
#include "quadrature.hpp"

template <class xType>
class location
{
public:
	typedef xType x_t;
	typedef typename x_t::Scalar scalar_t;

	/**
	 * \brief default constructor needed for elem_accelerator
	 */
	location()
	{
	}

	template <class elem_t>
	location(elem_t const &elem, quadrature_elem<typename elem_t::domain_t> const &q)
	{
		static_assert(std::is_same<x_t, typename elem_t::x_t>::value,
			"Element and descriptor location types must match");
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
	x_t m_x;
	scalar_t m_jacobian;
};


template <class xType>
class location_with_normal : public location<xType>
{
public:
	typedef xType x_t;

	/**
	 * \brief default constructor needed for elem_accelerator
	 */
	location_with_normal()
	{
	}

	template <class elem_t>
	location_with_normal(elem_t const &elem, quadrature_elem<typename elem_t::domain_t> const &q)
		: location<xType>(elem, q)
	{
		m_normal = elem.get_normal(q.get_xi());
		m_normal /= m_normal.norm();
	}

	x_t const &get_normal(void) const
	{
		return m_normal;
	}

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
	x_t m_normal;
};


#endif

