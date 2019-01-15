#ifndef NIHU_NEARLY_SINGULAR_COLLOCATIONAL_HPP_INCLUDED
#define NIHU_NEARLY_SINGULAR_COLLOCATIONAL_HPP_INCLUDED

#include "line_quad_store.hpp"
#include "location_normal.hpp"
#include "../core/element.hpp"
#include "../core/field.hpp"
#include "../core/inverse_mapping.hpp"
#include "../core/kernel.hpp"
#include "../core/nearly_singular_planar_constant_collocation_shortcut.hpp"
#include "../core/shapeset.hpp"
#include "../util/block_product.hpp"

namespace NiHu
{
template <class TrialField, class Kernel, 
	unsigned RadialOrder, unsigned TangentialOrder, 
	class Enable = void>
class nearly_singular_collocational;

template <class TrialField, class Kernel,
	unsigned RadialOrder, unsigned TangentialOrder>
class nearly_singular_collocational<TrialField, Kernel, RadialOrder, TangentialOrder,
	typename std::enable_if<
	element_traits::is_surface_element<typename TrialField::elem_t>::value
	>::type>
{
public:
	/** \brief quadrature orders stored as internal constants */
	enum {
		radial_order = RadialOrder,			/**< \brief radial quadr. order */
		tangential_order = TangentialOrder	/**< \brief tangential quadr order */
	};

	/** \brief the trial field type */
	typedef TrialField trial_field_t;
	/** \brief the element type */
	typedef typename trial_field_t::elem_t elem_t;

	/** \brief the original reference domain type */
	typedef typename elem_t::domain_t domain_t;
	/** \brief the reference coordinate vector type */
	typedef typename domain_t::xi_t xi_t;
	/** \brief the geometrical scalar type */
	typedef typename elem_t::scalar_t scalar_t;

	/** \brief the physical coordinate vector type */
	typedef typename elem_t::x_t x_t;

	/** \brief shape function set type */
	typedef typename trial_field_t::nset_t trial_nset_t;
	/** \brief the shape function vector type */
	typedef typename trial_nset_t::shape_t trial_n_shape_t;

	/** \brief the kernel type */
	typedef Kernel kernel_t;

	/** \brief the kernel's test input type */
	typedef typename kernel_traits<kernel_t>::test_input_t test_input_t;
	/** \brief the kernel's trial input type */
	typedef typename kernel_traits<kernel_t>::trial_input_t trial_input_t;

	/** \brief the kernel's weighted trial input type */
	typedef typename weighted_input<trial_input_t, elem_t>::type w_trial_input_t;

	/** \brief value type of the integral result */
	typedef typename semi_block_product_result_type<
		typename kernel_t::result_t, trial_n_shape_t
	>::type total_result_t;

	/** \brief constructor
	 * \param [in] elem the element
	 * \param [in] kernel the kernel
	 */
	nearly_singular_collocational(
		field_base<trial_field_t> const &trial_field,
		kernel_base<kernel_t> const &kernel)
		: m_elem(trial_field.get_elem())
		, m_kernel(kernel)
	{
	}

	template <class result_t>
	void integrate(result_t &&I, test_input_t const &tsi)
	{
		// the reference point
		x_t x0 = tsi.get_x();

		// perform inverse mapping to obtain image of reference point on element
		double tol = 1e-5;
		unsigned max_iter = 100;
		inverse_mapping<elem_t> im(m_elem);
		if (!im.eval(x0, tol, max_iter))
			throw std::runtime_error("Could not perform inverse mapping");
		auto res = im.get_result();
		m_xi0 = res.topRows(xi_t::RowsAtCompileTime);
		m_zeta = res(xi_t::RowsAtCompileTime);

		// compute linearized element
		elem_t lin_elem = m_elem.get_linearized_elem(m_xi0);

		// get jacobian at the reference image
		auto N0 = trial_nset_t::template eval_shape<0>(m_xi0);

		/** \todo check if plane_elem_helper_mid can be used here */
		// geometrical parameters (planar helpers)
		unsigned const N = domain_t::num_corners;
		for (unsigned n = 0; n < N; ++n)
		{
			xi_t c1 = domain_t::get_corner(n);			// corner
			xi_t c2 = domain_t::get_corner((n + 1) % N);	// next corner
			xi_t l = xi_t::Zero();
			l(0) = 1.0;
			if ((c2 - c1).norm() > 1e-3)
				l = (c2 - c1).normalized();					// side unit vector

			xi_t d1 = c1 - m_xi0;			// vector to corners
			xi_t d0 = d1 - l * d1.dot(l);		// perpendicular to side

			m_theta_lim[n] = std::atan2(d1(1), d1(0));	// corner angle
			m_theta0[n] = std::atan2(d0(1), d0(0));	// mid angle
			m_ref_distance[n] = d0.norm();			// distance to side
		}

		// iterate through triangles
		for (unsigned n = 0; n < N; ++n)
		{
			// get angular integration limits
			scalar_t t1 = m_theta_lim[n];
			scalar_t t2 = m_theta_lim[(n + 1) % N];

			/** todo this is for highly distorted elements */
			if (std::abs(t2 - t1) < 1e-3)
				continue;

			// we assume that the domain's corners are listed in positive order
			if (std::abs(t2 - t1) > M_PI)
				t2 += 2.0 * M_PI;

			// theta integration
			for (auto const &q_theta : line_quad_store<tangential_order>::quadrature)
			{
				// compute theta base point and weight
				scalar_t xx = q_theta.get_xi()(0);
				scalar_t theta = ((1.0 - xx) * t1 + (1.0 + xx) * t2) / 2.0;
				scalar_t w_theta = q_theta.get_w() * (t2 - t1) / 2.0;

				// reference domain's limit
				scalar_t rho_lim = m_ref_distance[n] / std::cos(theta - m_theta0[n]);

				// radial part of surface integration
				for (auto const &q_rho : line_quad_store<radial_order>::quadrature)
				{
					// quadrature base point and weight
					scalar_t rho = (1.0 + q_rho.get_xi()(0)) * rho_lim / 2.0;
					scalar_t w_rho = q_rho.get_w() * rho_lim / 2.0;

					// compute location in reference domain
					xi_t xi(rho*std::cos(theta), rho*std::sin(theta));
					xi += m_xi0;

					// evaluate weighted trial input
					w_trial_input_t tri(m_elem, xi);
					w_trial_input_t tri_lin(lin_elem, xi);

					// evaluate kernel
					typename kernel_t::result_t GJ = m_kernel(tsi, tri)
						* tri.get_jacobian();
					typename kernel_t::result_t GJ_lin = m_kernel(tsi, tri_lin)
						* tri_lin.get_jacobian();

					// get shape function
					auto N = trial_nset_t::template eval_shape<0>(xi);
					total_result_t F = semi_block_product(GJ, N);
					total_result_t F_lin = semi_block_product(GJ_lin, N0);

					I += w_theta * w_rho * rho * (F - F_lin);
				} // end of loop over radial nodes
			} // end of loop over tangential nodes
		} // end of loop over triangles

		typename kernel_t::result_t anal_res = 
			nearly_singular_planar_constant_collocation_shortcut<kernel_t, elem_t>::eval(
			tsi, lin_elem, m_kernel
		);

		I += semi_block_product(anal_res, N0);

	} // end of function integrate

private:
	elem_t const &m_elem;		/**< \brief the element reference */
	kernel_base<kernel_t> const &m_kernel;	/**< \brief the kernel reference */
	xi_t m_xi0;
	double m_zeta;

	double m_theta_lim[domain_t::num_corners];
	double m_theta0[domain_t::num_corners];
	double m_ref_distance[domain_t::num_corners];
};

} // end of namespace NiHu

#endif /* NIHU_NEARLY_SINGULAR_COLLOCATIONAL_HPP_INCLUDED */
