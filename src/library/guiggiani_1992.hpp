/** \file guiggiani_1992.hpp
 * \brief Guiggiani's method for hypersingular collocational integrals
 */

#ifndef GUIGGIANI_1992_HPP_INCLUDED
#define GUIGGIANI_1992_HPP_INCLUDED

#include <cmath>
#include "../util/crtp_base.hpp"
#include "../util/store_pattern.hpp"
#include "../core/field.hpp"
#include "../core/kernel.hpp"
#include "../core/gaussian_quadrature.hpp"
#include "../core/duffy_quadrature.hpp"
#include "location_normal.hpp"

/**
 * \brief acclerator class for storing duffy-type surface quadratures
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField, unsigned order>
class guiggiani_surface_accelerator
{
	// CRTP check
	static_assert(std::is_base_of<field_base<TestField>, TestField>::value,
		"The test field must be derived from field_base<TestField>");
	static_assert(std::is_base_of<field_base<TrialField>, TrialField>::value,
		"The trial field must be derived from field_base<TrialField>");
public:
	/** \brief template argument as nested type */
	typedef TestField test_field_t;
	/** \brief template argument as nested type */
	typedef TrialField trial_field_t;

	/** \brief test N-set type */
	typedef typename test_field_t::nset_t test_nset_t;

	/** \brief the duffy quadrature type */
	typedef duffy_quadrature<
		gauss_family_tag,
		typename trial_field_t::elem_t::lset_t
	> trial_duffy_t;

	/** \brief the stored quadrature type */
	typedef typename trial_duffy_t::ret_quad_type quadrature_t;

	/**
	 * \brief return the stored quadrature for a singular point
	 * \param [in] idx the singular point index
	 * \return quadrature singular at the specified point
	 */
	quadrature_t const &get_quadrature(unsigned idx) const
	{
		return m_quadratures[idx];
	}

	/**
	 * \brief constructor initialising the quadratures
	 */
	guiggiani_surface_accelerator(void)
	{
		for (unsigned idx = 0; idx < test_nset_t::num_nodes; ++idx)
			m_quadratures[idx] += trial_duffy_t::on_face(
			order, test_nset_t::corner_at(idx));
	}

private:
	/** \brief array of singular quadratures */
	quadrature_t m_quadratures[test_nset_t::num_nodes];
};


/**
 * \brief accelerator class for storing line quadratures
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 * \tparam order the quadrature order
 */
template <class TestField, class TrialField, unsigned order>
class guiggiani_line_accelerator
{
	// CRTP check
	static_assert(std::is_base_of<field_base<TestField>, TestField>::value,
	"The test field must be derived from field_base<TestField>");
	static_assert(std::is_base_of<field_base<TrialField>, TrialField>::value,
		"The trial field must be derived from field_base<TrialField>");
public:
	/** \brief template argument as nested type */
	typedef TestField test_field_t;
	/** \brief template argument as nested type */
	typedef TrialField trial_field_t;

	/** \brief test N-set type */
	typedef typename test_field_t::nset_t test_nset_t;
	/** \brief trial N-set type */
	typedef typename trial_field_t::nset_t trial_nset_t;

	/** \brief the reference domain type */
	typedef typename trial_field_t::elem_t::domain_t domain_t;
	/** \brief the reference location type */
	typedef typename domain_t::xi_t xi_t;

	/** \brief the stored quadrature type */
	typedef gaussian_quadrature<line_domain> quadrature_t;

	/**
	 * \brief return the stored quadrature for a singular point
	 * \param [in] idx the singular point index
	 * \return quadrature singular at the specified point
	 */
	quadrature_t const *get_quadrature(unsigned idx) const
	{
		return m_quadratures[idx];
	}

	/**
	 * \brief constructor initialising the quadratures
	 */
	guiggiani_line_accelerator(void)
	{
		quadrature_t source(order);

		for (unsigned idx = 0; idx < test_nset_t::num_nodes; ++idx)
		{
			Eigen::Matrix<double, 2, 1> coords;
			xi_t xi0 = test_nset_t::corner_at(idx);

			for (unsigned n = 0; n < domain_t::num_corners; ++n)
			{
				xi_t d1 = domain_t::get_corner(n)-xi0;
				xi_t d2 = domain_t::get_corner((n + 1) % domain_t::num_corners) - xi0;

				coords(0) = std::atan2(d1(1), d1(0));
				coords(1) = std::atan2(d2(1), d2(0));		// second angle limit
				if (coords(0) > coords(1))			// unwrap angles if needed
					coords(1) += 2.0 * M_PI;
				m_quadratures[idx][n] = source.transform<line_1_shape_set>(coords);
			}
		}
	}

private:
	/** \brief array of line quadratures */
	quadrature_t m_quadratures[test_nset_t::num_nodes][domain_t::num_corners];
};


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
	enum { quadrature_order = traits_t::order };

	typedef typename traits_t::test_field_t test_field_t;
	typedef typename traits_t::trial_field_t trial_field_t;
	typedef typename trial_field_t::elem_t elem_t;
	static_assert(std::is_same<elem_t, typename test_field_t::elem_t>::value,
		"Test and Trial element types must be the same for collocational Guiggiani integration");
	typedef typename trial_field_t::domain_t domain_t;
	typedef typename domain_t::xi_t xi_t;
	typedef typename elem_t::x_t x_t;
	typedef typename elem_t::scalar_t scalar_t;
	typedef typename trial_field_t::nset_t trial_nset_t;
	typedef typename test_field_t::nset_t test_nset_t;
	typedef typename trial_nset_t::shape_t trial_n_shape_t;

	typedef typename traits_t::kernel_t kernel_t;

	typedef typename kernel_traits<kernel_t>::test_input_t test_input_t;
	typedef typename kernel_traits<kernel_t>::trial_input_t trial_input_t;

	typedef typename merge<
		trial_input_t,
		typename build<normal_jacobian<typename trial_input_t::space_t> >::type
	>::type w_trial_input_t;

	typedef store<guiggiani_surface_accelerator<test_field_t, trial_field_t, quadrature_order> > surf_accel_store_t;
	typedef store<guiggiani_line_accelerator<test_field_t, trial_field_t, quadrature_order> > line_accel_store_t;

	typedef Eigen::Matrix<scalar_t, domain_t::dimension, domain_t::dimension> trans_t;

	/** \brief constructor
	 * \param [in] elem the element
	 * \param [in] kernel the kernel
	 */
	guiggiani_base(elem_t const &elem, kernel_t const &kernel)
		: m_elem(elem), m_kernel(kernel)
	{
	}

private:
	/** \brief store xi0-related quantities for acceleration
	 * \param [in] xi0 the singular point in the reference domain
	 */
	void compute_xi0(xi_t const &xi0)
	{
		m_xi0 = xi0;
		m_x0 = m_elem.get_x(m_xi0);
		typename elem_t::dx_t dx = m_elem.get_dx(m_xi0);
		double gamma = std::acos(dx.col(XI).dot(dx.col(ETA)) / dx.col(XI).norm() / dx.col(ETA).norm());

		m_T <<
			1.0, std::cos(gamma) * dx.col(ETA).norm() / dx.col(XI).norm(),
			0.0, std::sin(gamma) * dx.col(ETA).norm() / dx.col(XI).norm();

		m_Tinv = m_T.inverse();

		m_eta0 = m_T * m_xi0;
		
		m_J0_vector = m_elem.get_normal(m_xi0) / m_T.determinant();
		m_J0 = m_J0_vector.norm();
		m_n0 = m_J0_vector / m_J0;
		m_N0 = trial_nset_t::eval_shape(xi0);
	}

	/** \brief store theta-related quantities for acceleration
	 * \param [in] theta the angle in the reference domain
	 */
	void compute_theta(double theta)
	{
		xi_t eta(std::cos(theta), std::sin(theta));
		xi_t xi = m_Tinv * eta; // contains cos theta sin theta in xi system

		typename elem_t::dx_t dx = m_elem.get_dx(m_xi0);
		typename elem_t::ddx_t ddx = m_elem.get_ddx(m_xi0);

		m_A_vector = dx.col(XI) * xi(0) + dx.col(ETA) * xi(1);
		m_A = m_A_vector.norm();
		m_B_vector = ddx.col(XIXI) * xi(0)*xi(0) / 2.0 + ddx.col(XIETA) * xi(0)*xi(1) + ddx.col(ETAETA) * xi(1)*xi(1) / 2.0;

		m_J1_vector = (xi(0) * (ddx.col(XIXI).cross(dx.col(ETA)) + dx.col(XI).cross(ddx.col(XIETA))) +
			xi(1) * (ddx.col(ETAXI).cross(dx.col(ETA)) + dx.col(XI).cross(ddx.col(ETAETA)))) / m_T.determinant();

		typename trial_nset_t::dshape_t dN = trial_nset_t::eval_dshape(m_xi0);
		m_N1 = dN.col(XI)*xi(0) + dN.col(ETA)*xi(1);
	}

	/** \brief evaluate the surface integral
	 * \tparam result_t the result type
	 * \tparam quadrature_t the quadrature type
	 * \param [out] I the row where the integral is collected
	 */
	template <class result_t, class quadrature_t>
	void surface_integral(result_t I, quadrature_t const &quadrature)
	{
		// bind the kernel at the test input
		test_input_t test_input(m_elem, m_xi0);
		auto bound = m_kernel.bind(test_input);

		// iterate through the quadrature points
		for (auto it = quadrature.begin(); it != quadrature.end(); ++it)
		{
			// vector from the singular point in the reference domain
			xi_t const &xi = it->get_xi();				// tip
			xi_t eta = m_T * xi;
			auto deta = eta - m_eta0;
			double theta = std::atan2(deta(1), deta(0));	// angle
			double rho = deta.norm();					// length

			// compute theta-related quantities
			compute_theta(theta);
			derived().compute_Fm1Fm2();

			// evaluate the kernel
			w_trial_input_t trial_input(m_elem, xi);
			auto F = (
				bound(trial_input) * trial_input.get_jacobian() / m_T.determinant() *
				trial_nset_t::eval_shape(xi)
				).eval();

			// subtract the analytical singularity
			auto singular_part = (m_Fm2 / rho + m_Fm1) / (rho*rho);
			for (int j = 0; j < F.cols(); ++j)	// loop needed for scalar casting
				F(j) -= singular_part(j);

			// integrate the remaining regular part
			I += it->get_w() * m_T.determinant() * F;
		}
	}

	/** \brief compute line integrals
	 * \tparam result_t the result type
	 * \param [out] I the row where the integral is collected
	 */
	template <class result_t>
	void line_integrals(result_t I)
	{
		unsigned const N = domain_t::num_corners;

		gaussian_quadrature<line_domain> quad(quadrature_order);

		for (unsigned n = 0; n < N; ++n)
		{
			xi_t c1 = m_T * domain_t::get_corner(n);			// corner
			xi_t c2 = m_T * domain_t::get_corner((n + 1) % N);	// next corner
			xi_t l = (c2 - c1).normalized();					// side vector

			xi_t d2 = c2 - m_eta0;	// vectors to corners
			xi_t d1 = c1 - m_eta0;	// vectors to corners

			xi_t d0 = d2 - l*d2.dot(l);				// perpendicular to side

			double t1 = std::atan2(d1(1), d1(0));
			double t2 = std::atan2(d2(1), d2(0));
			if (t1 > t2)
				t2 += 2.0 * M_PI;
			double t0 = std::atan2(d0(1), d0(0));		// mid angle
			double d = d0.norm();					// distance to side

			// iterate through quadrature points
			for (auto it = quad.begin(); it != quad.end(); ++it)
			{
				double xx = it->get_xi()(0);
				double ww = it->get_w();

				double theta = (1.0 - xx) / 2.0*t1 + (1.0 + xx) / 2.0*t2;
				double w = ww * (t2 - t1) / 2.0;

				compute_theta(theta);
				derived().compute_Fm1Fm2();

				double rho_lim = d / std::cos(theta - t0);	// distance from xi0 to xi
				auto toadd = w * (m_Fm1 * std::log(rho_lim) - m_Fm2 / rho_lim);

				for (int j = 0; j < I.cols(); ++j)	// loop needed for scalar casting
					I(j) += toadd(j);
			}
		}
	}


public:
	template <class result_t>
	void integral(result_t &result)
	{
		for (unsigned idx = 0; idx < test_nset_t::num_nodes; ++idx)
		{
			compute_xi0(test_nset_t::corner_at(idx));
			surface_integral(result.row(idx), surf_accel_store_t::m_data.get_quadrature(idx));
			line_integrals(result.row(idx));
		}
	}

protected:
	elem_t const &m_elem;	/**< \brief the element reference */
	kernel_t const &m_kernel;	/**< \brief the kernel reference */

	xi_t m_xi0;				/**< \brief the source local coordinate */
	x_t m_x0;				/**< \brief the source point */

	trans_t m_T;
	trans_t m_Tinv;
	xi_t m_eta0;				/**< \brief the source local coordinate */

	x_t m_J0_vector;		/**< \brief the Jacobian vector at the source point */
	scalar_t m_J0;			/**< \brief the Jacobian at the source point */
	x_t m_n0;				/**< \brief the unit normal vector at the source point */
	trial_n_shape_t m_N0;	/**< \brief the shape function vector at the source point */

	x_t m_A_vector;			/**< \brief the location derivative vector */
	scalar_t m_A;			/**< \brief the magnitude of the location derivative */
	x_t m_B_vector;			/**< \brief the location second derivative vector */
	x_t m_J1_vector;		/**< \brief the linear part of the Jacobian vector */
	trial_n_shape_t m_N1;	/**< \brief the linear part of the shape function vector */

	trial_n_shape_t m_Fm2;	/**< \brief the -2-nd order Laurent coefficient */
	trial_n_shape_t m_Fm1;	/**< \brief the -1-st order Laurent coefficient */
};

// forward declaration
template <class TestField, class TrialField, class Kernel, unsigned order>
class guiggiani;

/** \brief traits of a guiggiani class */
template <class TestField, class TrialField, class Kernel, unsigned Order>
struct guiggiani_traits<guiggiani<TestField, TrialField, Kernel, Order> >
{
	typedef TestField test_field_t;
	typedef TrialField trial_field_t;
	typedef Kernel kernel_t;
	enum { order = Order };
};

#include "../library/laplace_kernel.hpp"

template <class TestField, class TrialField, unsigned order>
class guiggiani<TestField, TrialField, laplace_3d_HSP_kernel, order>
	: public guiggiani_base<guiggiani<TestField, TrialField, laplace_3d_HSP_kernel, order> >
{
	typedef guiggiani_base<guiggiani> base_t;
	typedef typename  base_t::kernel_t kernel_t;
	typedef typename base_t::elem_t elem_t;

public:
	guiggiani(elem_t const &elem, kernel_t const &kernel) : base_t(elem, kernel)
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

#include "../library/helmholtz_kernel.hpp"

template <class TestField, class TrialField, unsigned order>
class guiggiani<TestField, TrialField, helmholtz_3d_HSP_kernel<double>, order >
	: public guiggiani_base<guiggiani<TestField, TrialField, helmholtz_3d_HSP_kernel<double>, order > >
{
	typedef guiggiani_base<guiggiani> base_t;
	typedef typename  base_t::kernel_t kernel_t;
	typedef typename base_t::elem_t elem_t;

public:
	guiggiani(elem_t const &elem, kernel_t const &kernel) : base_t(elem, kernel)
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
