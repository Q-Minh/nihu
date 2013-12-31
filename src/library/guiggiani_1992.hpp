/** \file guiggiani_1992.hpp
 * \brief Guiggiani's method for hypersingular collocational integrals
 */

#ifndef GUIGGIANI_1992_HPP_INCLUDED
#define GUIGGIANI_1992_HPP_INCLUDED

#include <cmath>
#include "../core/field.hpp"

/** \brief traits of a guiggiani class
 * \tparam Derived the CRTP derived class
 */
template <class Derived>
struct guiggiani_traits;

/** \brief CRTP base class of all guiggiani classes
 * \tparam Derived the CRTP derived class
 */
template <class Derived>
class guiggiani_base
{
public:
	enum {XI = 0, ETA = 1, XIXI = 0, XIETA = 1, ETAXI = 1, ETAETA = 2};
	typedef guiggiani_traits<Derived> traits_t;

	typedef typename traits_t::field_t field_t;
	typedef typename field_t::elem_t elem_t;
	typedef typename elem_t::xi_t xi_t;
	typedef typename elem_t::x_t x_t;
	typedef typename elem_t::scalar_t scalar_t;
	typedef typename field_t::nset_t nset_t;
	typedef typename nset_t::shape_t n_shape_t;

	/** \brief constructor
	 * \param [in] elem the element
	 */
	guiggiani_base(elem_t const &elem) : m_elem(elem)
	{
	}

	/** \brief compute the coefficients at a local coordinate and an angle
	 * \param [in] xi0 the local coordinate
	 * \param [theta] theta the angle parameter
	 */
	void compute(xi_t const &xi0, double theta)
	{
		m_xi0 = xi0;

		auto c = std::cos(theta);
		auto s = std::sin(theta);

		auto dx = m_elem.get_dx(xi0);
		auto ddx = m_elem.get_ddx(xi0);

		m_A_vector = dx.col(XI) * c + dx.col(ETA) * s;
		m_A = m_A_vector.norm();
		m_B_vector = ddx.col(XIXI) * c*c / 2.0 + ddx.col(XIETA) * c*s + ddx.col(ETAETA) * s*s / 2.0;

		m_J0_vector = m_elem.get_normal(xi0);
		m_J1_vector = c * (ddx.col(XIXI).cross(dx.col(ETA)) + dx.col(XI).cross(ddx.col(XIETA))) +
			s * (ddx.col(ETAXI).cross(dx.col(ETA)) + dx.col(XI).cross(ddx.col(ETAETA)));
		m_J0 = m_J0_vector.norm();
		m_n0 = m_J0_vector / m_J0;

		m_N0 = nset_t::eval_shape(xi0);
		auto dN = nset_t::eval_dshape(xi0);
		m_N1 = dN.col(XI)*c + dN.col(ETA)*s;
	}

protected:
	elem_t const &m_elem;	/**< \brief the element reference */
	xi_t m_xi0;				/**< \brief the local coordinate */
	x_t m_A_vector;			/**< \brief the location derivative vector */
	scalar_t m_A;			/**< \brief the magnitude of the location derivative */
	x_t m_B_vector;			/**< \brief the location second derivative vector */
	x_t m_J0_vector;		/**< \brief the constant Jacobian vector */
	scalar_t m_J0;			/**< \brief the magnitude of the constant Jacobian */
	x_t m_n0;				/**< \brief the unit normal vector */
	x_t m_J1_vector;		/**< \brief the linear part of the Jacobian vector */
	n_shape_t m_N0;			/**< \brief the constant part of the shape function vector */
	n_shape_t m_N1;			/**< \brief the linear part of the shape function vector */
	n_shape_t m_Fm2;		/**< \brief the -2-nd order Laurent coefficient */
	n_shape_t m_Fm1;		/**< \brief the -1-st order Laurent coefficient */
};

// forward declaration
template <class Field, class Kernel>
class guiggiani;

/** \brief traits of a guiggiani class */
template <class Field, class Kernel>
struct guiggiani_traits<guiggiani<Field, Kernel> >
{
	typedef Field field_t;
	typedef Kernel kernel_t;
};

#include "../library/laplace_kernel.hpp"

template <class Field>
class guiggiani<Field, laplace_3d_HSP_kernel> : public guiggiani_base<guiggiani<Field, laplace_3d_HSP_kernel> >
{
public:
	typedef guiggiani_base<guiggiani> base_t;
	typedef typename base_t::elem_t elem_t;

	guiggiani(elem_t const &elem) : base_t(elem)
	{
	}

	void compute_Fm1Fm2(void)
	{
		auto A2 = this->m_A*this->m_A;
		auto A3 = A2 * this->m_A;
		this->m_Fm2 = this->m_J0 * this->m_N0 / (4.0 * M_PI * A3);
		this->m_Fm1 = (this->m_J0*this->m_N1 + this->m_J1_vector.dot(this->m_n0)*this->m_N0
			- 3.0*this->m_N0*this->m_J0* this->m_A_vector.dot(this->m_B_vector) / A2) / (4.0*M_PI * A3);
	}
};

#endif // GUIGGIANI_1992_HPP_INCLUDED
