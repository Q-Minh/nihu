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

	/** \brief the traits type */
	typedef guiggiani_traits<Derived> traits_t;

	enum { XI = 0, ETA = 1, XIXI = 0, XIETA = 1, ETAXI = 1, ETAETA = 2 };

	/** \brief the quadrature order */
	enum { quadrature_order = traits_t::order };

	/** \brief the test field type */
	typedef typename traits_t::test_field_t test_field_t;
	/** \brief the trial field type */
	typedef typename traits_t::trial_field_t trial_field_t;
	/** \brief the element type */
	typedef typename trial_field_t::elem_t elem_t;

	// compile time check that the test and trial elements are of the same type
	static_assert(std::is_same<elem_t, typename test_field_t::elem_t>::value,
		"Test and Trial element types must be the same for collocational Guiggiani integration");

	/** \brief the reference domain type */
	typedef typename elem_t::domain_t domain_t;
	/** \brief the reference coordinate vector type */
	typedef typename domain_t::xi_t xi_t;
	/** \brief the physical coordinate vector type */
	typedef typename elem_t::x_t x_t;
	/** \brief the scalar type */
	typedef typename elem_t::scalar_t scalar_t;

	/** \brief the trial N-set type */
	typedef typename trial_field_t::nset_t trial_nset_t;
	/** \brief the test N-set type */
	typedef typename test_field_t::nset_t test_nset_t;
	/** \brief the trial shape function vector type */
	typedef typename trial_nset_t::shape_t trial_n_shape_t;

	/** \brief the kernel type */
	typedef typename traits_t::kernel_t kernel_t;

	/** \brief the kernel's test input type */
	typedef typename kernel_traits<kernel_t>::test_input_t test_input_t;
	/** \brief the kernel's trial input type */
	typedef typename kernel_traits<kernel_t>::trial_input_t trial_input_t;

	/** \brief the kernel's weighted trial input type */
	typedef typename merge<
		trial_input_t,
		typename build<normal_jacobian<typename trial_input_t::space_t> >::type
	>::type w_trial_input_t;

	/** \brief the Rong's transformation matrix type */
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
		scalar_t gamma = std::acos(dx.col(XI).dot(dx.col(ETA)) / dx.col(XI).norm() / dx.col(ETA).norm());

		m_T <<
			1.0, std::cos(gamma) * dx.col(ETA).norm() / dx.col(XI).norm(),
			0.0, std::sin(gamma) * dx.col(ETA).norm() / dx.col(XI).norm();

		m_Tinv = m_T.inverse();

		m_eta0 = m_T * m_xi0;

		m_J0_vector = m_elem.get_normal(m_xi0) / m_T.determinant();
		m_J0 = m_J0_vector.norm();
		m_n0 = m_J0_vector / m_J0;
		m_N0 = trial_nset_t::eval_shape(xi0);

		// geometrical parameters (triangle helpers)
		unsigned const N = domain_t::num_corners;
		for (unsigned n = 0; n < N; ++n)
		{
			xi_t c1 = m_T * domain_t::get_corner(n);			// corner
			xi_t c2 = m_T * domain_t::get_corner((n + 1) % N);	// next corner
			xi_t l = (c2 - c1).normalized();					// side vector

			xi_t d1 = c1 - m_eta0;			// vector to corners
			xi_t d0 = d1 - l*d1.dot(l);		// perpendicular to side

			m_theta_lim[n] = std::atan2(d1(1), d1(0));	// corner angle
			m_theta0[n] = std::atan2(d0(1), d0(0));		// mid angle
			m_ref_distance[n] = d0.norm();				// distance to side
		}

		// print_debug_xi0();
	}

	void print_debug_xi0(void)
	{
		std::cout << "PrintDebug xi0" << std::endl;
		std::cout << "xi0:\t" << m_xi0.transpose() << std::endl;
		std::cout << "x0:\t" << m_x0.transpose() << std::endl;
		std::cout << "T:\t" << m_T << std::endl;
		std::cout << "T(-1):\t" << m_Tinv << std::endl;
		std::cout << "eta0:\t" << m_eta0.transpose() << std::endl;
		std::cout << "J0vec:\t" << m_J0_vector.transpose() << std::endl;
		std::cout << "J0:\t" << m_J0 << std::endl;
	}

	/** \brief store theta-related quantities for acceleration
	 * \param [in] theta the angle in the reference domain
	 */
	void compute_theta(scalar_t theta)
	{
		// contains cos theta sin theta in xi system
		xi_t xi = m_Tinv.col(0) * std::cos(theta) + m_Tinv.col(1) * std::sin(theta);

		typename elem_t::dx_t dx = m_elem.get_dx(m_xi0);
		typename elem_t::ddx_t ddx = m_elem.get_ddx(m_xi0);

		m_A_vector = dx.col(XI) * xi(0) + dx.col(ETA) * xi(1);
		m_A = m_A_vector.norm();
		m_B_vector = ddx.col(XIXI) * xi(0)*xi(0) / 2.0 + ddx.col(XIETA) * xi(0)*xi(1) + ddx.col(ETAETA) * xi(1)*xi(1) / 2.0;

		m_J1_vector = (xi(0) * (ddx.col(XIXI).cross(dx.col(ETA)) + dx.col(XI).cross(ddx.col(XIETA))) +
			xi(1) * (ddx.col(ETAXI).cross(dx.col(ETA)) + dx.col(XI).cross(ddx.col(ETAETA)))) / m_T.determinant();

		auto dN = trial_nset_t::eval_dshape(m_xi0);
		m_N1 = dN.col(XI)*xi(0) + dN.col(ETA)*xi(1);
	}

	void print_debug_theta(scalar_t theta)
	{
		std::cout << "\nPrintDebug theta\n" << std::endl;
		std::cout << "theta:\t" << theta << std::endl;
		std::cout << "Avec:\t" << m_A_vector.transpose() << std::endl;
		std::cout << "A:\t" << m_A << std::endl;
		std::cout << "Bvec:\t" << m_B_vector.transpose() << std::endl;
		std::cout << "J1vec:\t" << m_J1_vector.transpose() << std::endl;
		std::cout << "F1:\t" << m_Fm1 << std::endl;
		std::cout << "F2:\t" << m_Fm2 << std::endl;
	}

	/** \brief evaluate the surface integral
	 * \tparam result_t the result type
	 * \param [out] I the row where the integral is collected
	 */
	template <class result_t>
	void integrate(result_t I)
	{
		// bind the kernel at the test input
		test_input_t test_input(m_elem, m_xi0);
		auto bound = m_kernel.bind(test_input);

		gaussian_quadrature<line_domain> quad(quadrature_order);

		// iterate through triangles
		unsigned const N = domain_t::num_corners;
		for (unsigned n = 0; n < N; ++n)
		{
			scalar_t t1 = m_theta_lim[n];
			scalar_t t2 = m_theta_lim[(n + 1) % N];
			if (t2 < t1)
				t2 += 2.0 * M_PI;

			// theta integration
			for (auto it_theta = quad.begin(); it_theta != quad.end(); ++it_theta)
			{
				// compute theta base point and weight
				scalar_t xx = it_theta->get_xi()(0);
				scalar_t theta = (1.0 - xx) / 2.0 * t1 + (1.0 + xx) / 2.0 * t2;
				scalar_t w_theta = it_theta->get_w() * (t2 - t1) / 2.0;

				// compute theta-related members
				compute_theta(theta);
				derived().compute_Fm1Fm2();

				// reference domain's limit
				scalar_t rho_lim = m_ref_distance[n] / std::cos(theta - m_theta0[n]);	// distance from xi0 to xi

				auto toadd = w_theta * (m_Fm1 * std::log(rho_lim) - m_Fm2 / rho_lim);
				for (int j = 0; j < I.cols(); ++j)	// loop needed for scalar casting
					I(j) += toadd(j);

				// rho integration
				for (auto it_rho = quad.begin(); it_rho != quad.end(); ++it_rho)
				{
					// quadrature base point and weight
					scalar_t rho = (1.0 + it_rho->get_xi()(0)) * rho_lim / 2.0;
					scalar_t w_rho = it_rho->get_w() * rho_lim / 2.0;

					// compute location in original reference domain
					xi_t eta(rho*std::cos(theta), rho*std::sin(theta));
					eta += m_eta0;
					xi_t xi = m_Tinv * eta;

					// evaluate kernel multiplied by rho and the jacobian and shape vector
					w_trial_input_t trial_input(m_elem, xi);
					auto F = (
						bound(trial_input) *
						trial_input.get_jacobian() / m_T.determinant() *
						rho *
						trial_nset_t::eval_shape(xi)
						).eval();

					// subtract the analytical singularity
					auto singular_part = (m_Fm2 / rho + m_Fm1) / rho;
					for (int j = 0; j < F.cols(); ++j)	// loop needed for scalar casting
						F(j) -= singular_part(j);

					// integration
					I += w_theta * w_rho * F;
				} // end of rho loop
			} // end of theta loop
		} // end of loop on triangle sides
	} // end of function

public:
	template <class result_t>
	void integral(result_t &result)
	{
		for (unsigned idx = 0; idx < test_nset_t::num_nodes; ++idx)
		{
			compute_xi0(test_nset_t::corner_at(idx));
			integrate(result.row(idx));
		}
	}

protected:
	elem_t const &m_elem;	/**< \brief the element reference */
	kernel_t const &m_kernel;	/**< \brief the kernel reference */

	xi_t m_xi0;				/**< \brief the source local coordinate */
	x_t m_x0;				/**< \brief the source point */

	trans_t m_T;			/**< \brief transformation to distorted reference domain */
	trans_t m_Tinv;			/**< \brief inverse of the transformation */
	xi_t m_eta0;			/**< \brief the source local coordinate in the distorted reference domain */

	scalar_t m_theta_lim[domain_t::num_corners];
	scalar_t m_theta0[domain_t::num_corners];
	scalar_t m_ref_distance[domain_t::num_corners];

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

template <class TestField, class TrialField, class WaveNumber, unsigned order>
class guiggiani<TestField, TrialField, helmholtz_3d_HSP_kernel<WaveNumber>, order >
	: public guiggiani_base<guiggiani<TestField, TrialField, helmholtz_3d_HSP_kernel<WaveNumber>, order > >
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
