/** \file guiggiani_1992.hpp
 * \brief Guiggiani's method for hypersingular collocational integrals
 */

#ifndef GUIGGIANI_1992_HPP_INCLUDED
#define GUIGGIANI_1992_HPP_INCLUDED

#include <cmath>
#include "../core/shapeset.hpp"
#include "laplace_kernel.hpp"

/** \brief field-related quantities of the Guiggiani-method
 * \tparam Field the field type
 */
template <class Field>
class guiggiani_field
{
	/** \brief template argument as nested type */
	typedef Field field_t;

	/** \brief the element type */
	typedef typename field_t::elem_t elem_t;
	/** \brief the element L-set */
	typedef typename elem_t::lset_t lset_t;
	/** \brief the element N-set */
	typedef typename field_t::nset_t nset_t;

	/** \brief the intrinsic coordinate vctor type */
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

	static typename elem_t::x_t J0_vector(elem_t const &elem, xi_t const &xi0)
	{
		return elem.get_normal(xi0);
	}

	static typename elem_t::x_t J1_vector(elem_t const &elem, xi_t const &xi0, double theta)
	{
		auto dx = elem.get_dx(xi0);
		auto ddx = elem.get_ddx(xi0);

		return std::cos(theta) * (ddx.col(0).cross(dx.col(1)) + dx.col(0).cross(ddx.col(1))) +
			std::sin(theta) * (ddx.col(1).cross(dx.col(1)) + dx.col(0).cross(ddx.col(2)));
	}

	/** \brief evaluate static part of the shape functions
	 * \param [in] xi0 the reference coordinate
	 * \return the static part of the shape functions
	 */
	static typename nset_t::shape_t N0(xi_t const &xi0)
	{
		return nset_t::eval_shape(xi0);
	}

	/** \brief evaluate linear part of the shape functions
	 * \param [in] xi0 the reference coordinate
	 * \param [in] theta the angle parameter
	 * \return the linear part of the shape functions
	 */
	static typename nset_t::shape_t N1(xi_t const &xi0, double theta)
	{
		auto dN = nset_t::eval_dshape(xi0);
		return dN.col(0)*cos(theta) + dN.col(1)*sin(theta);
	}
};



template <class Kernel, class Field>
class guiggiani_hypersingular;

template <class Field>
class guiggiani_hypersingular<laplace_3d_HSP_kernel, Field>
{
	typedef Field field_t;

	typedef typename field_t::nset_t nset_t;
	typedef typename nset_t::shape_t n_shape_t;
	/** \brief the element type */
	typedef typename field_t::elem_t elem_t;
	/** \brief the intrinsic coordinate vctor type */
	typedef typename elem_t::xi_t xi_t;

	typedef guiggiani_field<field_t> guifield_t;

public:
	static n_shape_t Fm2(elem_t const &elem, xi_t const &xi0, double theta)
	{
		auto A = guifield_t::A_vector(elem, xi0, theta).norm();
		auto J0 = guifield_t::J0_vector(elem, xi0).norm();
		return J0 * guifield_t::N0(xi0) / (4.0 * M_PI * A*A*A);
	}

	static n_shape_t Fm1(elem_t const &elem, xi_t const &xi0, double theta)
	{
		auto Avec = guifield_t::A_vector(elem, xi0, theta);
		auto Bvec = guifield_t::B_vector(elem, xi0, theta);
		auto A = Avec.norm();
		Avec /= A;
		Bvec /= A;
		auto J0vec = guifield_t::J0_vector(elem, xi0);
		auto J0 = J0vec.norm();
		auto nx0 = J0vec / J0;
		auto J1 = guifield_t::J1_vector(elem, xi0, theta).dot(nx0);
		auto N0 = guifield_t::N0(xi0);
		auto N1 = guifield_t::N1(xi0, theta);

		return (J0*N1 + J1*N0 - 3.0*N0*J0*Avec.dot(Bvec)) / (4.0*M_PI * A*A*A);
	}

private:
};

#endif // GUIGGIANI_1992_HPP_INCLUDED
