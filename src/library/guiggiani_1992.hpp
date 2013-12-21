/** \file guiggiani_1992.hpp
 * \brief Guiggiani's method for hypersingular collocational integrals
 */

#ifndef GUIGGIANI_1992_HPP_INCLUDED
#define GUIGGIANI_1992_HPP_INCLUDED

#include <cmath>
#include "../core/shapeset.hpp"

template <class Kernel, class Field>
class guiggiani_hypersingular_collocation
{
	/** \brief template argument as nested type */
	typedef Field field_t;

	typedef typename field_t::elem_t elem_t;
	typedef typename elem_t::lset_t lset_t;
	typedef typename field_t::nset_t nset_t;

	typedef typename lset_t::shape_t l_shape_t;
	typedef typename nset_t::shape_t n_shape_t;

	typedef typename elem_t::xi_t xi_t;

public:

	/** \brief compute location derivative with respect to rho
	 * \param [in] xi0 the reference location at the collocation point
	 * \return location derivatives
	 */
	static typename elem_t::x_t A_vector(elem_t const &elem, xi_t const &xi0, double theta)
	{
		auto dx = elem.get_dx(xi0);
		return dx.col(0) * std::cos(theta) + dx.col(1) * std::sin(theta);
	}

	/** \brief compute location second derivative with respect to rho
	 * \param [in] xi0 the reference location at the collocation point
	 * \return location second derivative vector
	 */
	static typename elem_t::x_t B_vector(elem_t const &elem, xi_t const &xi0, double theta)
	{
		auto ddx = elem.get_ddx(xi0);
		auto c = std::cos(theta);
		auto s = std::sin(theta);
		return ddx.col(0) * c*c/2.0 + ddx.col(1) * c*s + ddx.col(2) * s*s/2.0;
	}

	/** \brief evaluate static part of the shape functions
	 * \param [in] xi0 the reference coordinate
	 * \return the static part of the shape functions
	 */
	static n_shape_t N0(xi_t const &xi0)
	{
		return nset_t::eval_shape(xi0);
	}

	/** \brief evaluate linear part of the shape functions
	 * \param [in] xi0 the reference coordinate
	 * \param [in] theta the angle parameter
	 * \return the linear part of the shape functions
	 */
	static n_shape_t N1(xi_t const &xi0, double theta)
	{
		auto dN = nset_t::eval_dshape(xi0);
		return dN.col(0)*cos(theta) + dN.col(1)*sin(theta);
	}
};

#endif // GUIGGIANI_1992_HPP_INCLUDED
