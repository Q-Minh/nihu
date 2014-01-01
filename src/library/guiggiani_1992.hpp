/** \file guiggiani_1992.hpp
 * \brief Guiggiani's method for hypersingular collocational integrals
 */

#ifndef GUIGGIANI_1992_HPP_INCLUDED
#define GUIGGIANI_1992_HPP_INCLUDED

#include <cmath>
#include "../util/crtp_base.hpp"
#include "../core/field.hpp"
#include "../core/gaussian_quadrature.hpp"


/** \brief store-wrapper of a statically stored line quadrature */
template <unsigned order>
struct line_quad_store
{
	/** \brief the stored static quadrature member */
	static gaussian_quadrature<line_domain> const quadrature;
};

/** \brief definition of the statically stored quadrature member */
template <unsigned order>
gaussian_quadrature<line_domain> const line_quad_store<order>::quadrature(order);


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
	NIHU_CRTP_HELPERS

	typedef guiggiani_traits<Derived> traits_t;

	enum { XI = 0, ETA = 1, XIXI = 0, XIETA = 1, ETAXI = 1, ETAETA = 2 };
	enum { quadrature_order = 7 };
	typedef line_quad_store<quadrature_order> line_quadr_t;

	typedef typename traits_t::field_t field_t;
	typedef typename field_t::elem_t elem_t;
	typedef typename field_t::domain_t domain_t;
	typedef typename domain_t::xi_t xi_t;
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

		m_J0_vector = m_elem.get_normal(xi0);
		m_J0 = m_J0_vector.norm();
		m_n0 = m_J0_vector / m_J0;
		m_N0 = nset_t::eval_shape(xi0);

		m_A_vector = dx.col(XI) * c + dx.col(ETA) * s;
		m_A = m_A_vector.norm();
		m_B_vector = ddx.col(XIXI) * c*c / 2.0 + ddx.col(XIETA) * c*s + ddx.col(ETAETA) * s*s / 2.0;

		m_J1_vector = c * (ddx.col(XIXI).cross(dx.col(ETA)) + dx.col(XI).cross(ddx.col(XIETA))) +
			s * (ddx.col(ETAXI).cross(dx.col(ETA)) + dx.col(XI).cross(ddx.col(ETAETA)));

		auto dN = nset_t::eval_dshape(xi0);
		m_N1 = dN.col(XI)*c + dN.col(ETA)*s;
	}

	void line_integrals(xi_t const &xi0, n_shape_t &I1, n_shape_t &I2)
	{
		unsigned const N = domain_t::num_corners;
		for (unsigned n = 0; n < N; ++n)
		{
			auto const &c1 = domain_t::get_corner(n);
			auto const &c2 = domain_t::get_corner((n + 1) % N);
			auto l = (c2 - c1).normalized();

			auto d1 = c1 - xi0, d2 = c2 - xi0;
			auto d0 = d2 - l*d2.dot(l);
			double t1 = std::atan2(d1(1), d1(0)), t2 = std::atan2(d2(1), d2(0));
			if (t1 > t2)
				t2 += 2.0 * M_PI;
			double t0 = std::atan2(d0(1), d0(0));
			double d = d0.norm();

			for (auto it = line_quadr_t::quadrature.begin(); it != line_quadr_t::quadrature.end(); ++it)
			{
				auto x = it->get_xi()(0);
				auto w = it->get_w() * (t2 - t1) / 2.0;
				auto theta = (t1 * (1 - x) + t2 * (1 + x)) / 2.0;

				compute(xi0, theta);
				derived().compute_Fm1Fm2();

				auto beta = 1.0 / m_A;
				auto gamma = -m_A_vector.dot(m_B_vector) / (m_A*m_A*m_A*m_A);

				double rho_lim = d / std::cos(theta-t0);

				I1 += w * m_Fm1 * std::log(std::abs(rho_lim / beta));
				I2 += w * m_Fm2 * (gamma / beta / beta + 1.0 / rho_lim);
			}
		}
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
