/** \file guiggiani_1992.hpp
 * \brief Guiggiani's method for hypersingular collocational integrals
 */

#ifndef GUIGGIANI_1992_HPP_INCLUDED
#define GUIGGIANI_1992_HPP_INCLUDED

#include <cmath>
#include "../util/crtp_base.hpp"
#include "../core/field.hpp"
#include "../core/kernel.hpp"
#include "../core/gaussian_quadrature.hpp"
#include "location_normal.hpp"

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


/** \brief store-wrapper of a statically stored surface quadrature */
template <class domain, unsigned order>
struct surface_quad_store
{
	/** \brief the stored static quadrature member */
	static gaussian_quadrature<domain> const quadrature;
};

/** \brief definition of the statically stored quadrature member */
template <class domain, unsigned order>
gaussian_quadrature<domain> const surface_quad_store<domain, order>::quadrature(order);


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

	typedef typename traits_t::field_t field_t;
	typedef typename field_t::elem_t elem_t;
	typedef typename field_t::domain_t domain_t;
	typedef typename domain_t::xi_t xi_t;
	typedef typename elem_t::x_t x_t;
	typedef typename elem_t::scalar_t scalar_t;
	typedef typename field_t::nset_t nset_t;
	typedef typename nset_t::shape_t n_shape_t;

	typedef typename traits_t::kernel_t kernel_t;

	typedef typename kernel_traits<kernel_t>::test_input_t test_input_t;
	typedef typename kernel_traits<kernel_t>::trial_input_t trial_input_t;

	typedef test_input_t w_test_input_t;
	typedef typename merge<
		trial_input_t,
		typename build<normal_jacobian<typename trial_input_t::space_t> >::type
	>::type w_trial_input_t;

	typedef line_quad_store<quadrature_order> line_quadr_t;
	typedef surface_quad_store<domain_t, quadrature_order> surf_quadr_t;

	/** \brief constructor
	 * \param [in] elem the element
	 */
	guiggiani_base(elem_t const &elem) : m_elem(elem)
	{
	}

private:
	void compute_xi0(xi_t const &xi0)
	{
		m_xi0 = xi0;
		m_x0 = m_elem.get_x(m_xi0);
		m_J0_vector = m_elem.get_normal(m_xi0);
		m_J0 = m_J0_vector.norm();
		m_n0 = m_J0_vector / m_J0;

		m_N0 = nset_t::eval_shape(xi0);
	}

	void compute_theta(double theta)
	{
		auto c = std::cos(theta);
		auto s = std::sin(theta);

		typename elem_t::dx_t dx = m_elem.get_dx(m_xi0);
		typename elem_t::ddx_t ddx = m_elem.get_ddx(m_xi0);

		m_A_vector = dx.col(XI) * c + dx.col(ETA) * s;
		m_A = m_A_vector.norm();
		m_B_vector = ddx.col(XIXI) * c*c / 2.0 + ddx.col(XIETA) * c*s + ddx.col(ETAETA) * s*s / 2.0;

		m_J1_vector = c * (ddx.col(XIXI).cross(dx.col(ETA)) + dx.col(XI).cross(ddx.col(XIETA))) +
			s * (ddx.col(ETAXI).cross(dx.col(ETA)) + dx.col(XI).cross(ddx.col(ETAETA)));

		typename nset_t::dshape_t dN = nset_t::eval_dshape(m_xi0);
		m_N1 = dN.col(XI)*c + dN.col(ETA)*s;
	}

	void surface_integral(n_shape_t &I)
	{
		for (auto it = surf_quadr_t::quadrature.begin(); it != surf_quadr_t::quadrature.end(); ++it)
		{
			xi_t const &xi = it->get_xi();
			xi_t dxi = xi - m_xi0;
			double theta = std::atan2(dxi(1), dxi(0));
			double rho = dxi.norm();

			compute_theta(theta);
			derived().compute_Fm1Fm2();
			
			auto F = derived().compute_kernel(xi);
			I += it->get_w() * (F - (m_Fm2 / rho + m_Fm1) / rho / rho);
		}
	}

	void line_integrals(n_shape_t &I)
	{
		unsigned const N = domain_t::num_corners;
		for (unsigned n = 0; n < N; ++n)
		{
			xi_t const &c1 = domain_t::get_corner(n);
			xi_t const &c2 = domain_t::get_corner((n + 1) % N);
			xi_t l = (c2 - c1).normalized();

			xi_t d1 = c1 - m_xi0, d2 = c2 - m_xi0;
			xi_t d0 = d2 - l*d2.dot(l);
			double t1 = std::atan2(d1(1), d1(0)), t2 = std::atan2(d2(1), d2(0));
			if (t1 > t2)
				t2 += 2.0 * M_PI;
			double t0 = std::atan2(d0(1), d0(0));
			double d = d0.norm();

			for (auto it = line_quadr_t::quadrature.begin(); it != line_quadr_t::quadrature.end(); ++it)
			{
				double x = it->get_xi()(0);
				double w = it->get_w() * (t2 - t1) / 2.0;
				double theta = (t1 * (1 - x) + t2 * (1 + x)) / 2.0;

				compute_theta(theta);
				derived().compute_Fm1Fm2();

				double rho_lim = d / std::cos(theta - t0);

				I += w * (m_Fm1 * std::log(rho_lim) - m_Fm2 / rho_lim);
			}
		}
	}

public:
	void integral(xi_t const &xi0, n_shape_t &result)
	{
		compute_xi0(xi0);
		surface_integral(result);
		line_integrals(result);
	}

protected:
	elem_t const &m_elem;	/**< \brief the element reference */

	xi_t m_xi0;				/**< \brief the source local coordinate */
	x_t m_x0;				/**< \brief the source point */
	x_t m_J0_vector;		/**< \brief the Jacobian vector at the source point */
	scalar_t m_J0;			/**< \brief the Jacobian at the source point */
	x_t m_n0;				/**< \brief the unit normal vector at the source point */
	n_shape_t m_N0;			/**< \brief the shape function vector at the source point */

	x_t m_A_vector;			/**< \brief the location derivative vector */
	scalar_t m_A;			/**< \brief the magnitude of the location derivative */
	x_t m_B_vector;			/**< \brief the location second derivative vector */
	x_t m_J1_vector;		/**< \brief the linear part of the Jacobian vector */
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
	typedef typename base_t::nset_t nset_t;
	typedef typename base_t::n_shape_t n_shape_t;

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

	n_shape_t compute_kernel(typename elem_t::xi_t const &xi)
	{
		auto y = this->m_elem.get_x(xi);
		auto Jny = this->m_elem.get_normal(xi);

		auto rvec = y - this->m_x0;
		auto r = rvec.norm();
		auto rdnx = -rvec.dot(this->m_n0) / r;
		auto rdny = rvec.dot(Jny) / r;

		auto N = nset_t::eval_shape(xi);
		return N / (4.0*M_PI * r*r*r) * (this->m_n0.dot(Jny) + 3.0 * rdnx * rdny);
	}
};


template <class Field>
class guiggiani<Field, laplace_2d_HSP_kernel> : public guiggiani_base<guiggiani<Field, laplace_2d_HSP_kernel> >
{
public:
	typedef guiggiani_base<guiggiani> base_t;
	typedef typename base_t::elem_t elem_t;
	typedef typename base_t::nset_t nset_t;
	typedef typename base_t::n_shape_t n_shape_t;

	guiggiani(elem_t const &elem) : base_t(elem)
	{
	}

	void compute_Fm1Fm2(void)
	{
		auto A2 = this->m_A*this->m_A;
		auto A3 = A2 * this->m_A;
		this->m_Fm2 = this->m_J0 * this->m_N0 / (2.0 * M_PI * A2);
		this->m_Fm1 = (this->m_J0*this->m_N1 + this->m_J1_vector.dot(this->m_n0)*this->m_N0
			- 2.0*this->m_N0*this->m_J0* this->m_A_vector.dot(this->m_B_vector) / A2) / (2.0*M_PI * A2);
	}

	n_shape_t compute_kernel(typename elem_t::xi_t const &xi)
	{
		auto y = this->m_elem.get_x(xi);
		auto Jny = this->m_elem.get_normal(xi);

		auto rvec = y - this->m_x0;
		auto r = rvec.norm();
		auto rdnx = -rvec.dot(this->m_n0) / r;
		auto rdny = rvec.dot(Jny) / r;

		auto N = nset_t::eval_shape(xi);
		return N / (2.0*M_PI * r*r) * (this->m_n0.dot(Jny) + 2.0 * rdnx * rdny);
	}
};


#endif // GUIGGIANI_1992_HPP_INCLUDED
