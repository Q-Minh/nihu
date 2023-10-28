#ifndef NIHU_NEARLY_SINGULAR_COLLOCATIONAL_TELLES_HPP_INCLUDED
#define NIHU_NEARLY_SINGULAR_COLLOCATIONAL_TELLES_HPP_INCLUDED

#include "line_quad_store.hpp"
#include "location_normal.hpp"
#include "weighted_input.hpp"
#include "telles_1987.hpp"
#include "../core/element.hpp"
#include "../core/field.hpp"
#include "../core/inverse_mapping.hpp"
#include "../core/kernel.hpp"
#include "../core/nearly_singular_planar_constant_collocation_shortcut.hpp"
#include "../core/shapeset.hpp"
#include "../util/block_product.hpp"

namespace NiHu
{
template <class TrialField, class Kernel, unsigned Order, class Enable = void>
class nearly_singular_collocational_telles;

template <class TrialField, class Kernel, unsigned Order>
class nearly_singular_collocational_telles<TrialField, Kernel, Order,
	typename std::enable_if<
	element_traits::is_surface_element<typename TrialField::elem_t>::value
	>::type>
{
public:
	/** \brief quadrature orders stored as internal constants */
	enum {
		order = Order,			/**< \brief quadr. order */
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

	typedef gaussian_quadrature<domain_t> quad_t;

	/** \brief constructor
	 * \param [in] trial_field the element
	 * \param [in] kernel the kernel
	 */
	nearly_singular_collocational_telles(
		field_base<trial_field_t> const &trial_field,
		kernel_base<kernel_t> const &kernel)
		: m_elem(trial_field.get_elem())
		, m_kernel(kernel)
		, m_quad_orig(Order)
	{
	}

	template <class result_t>
	void integrate(result_t &&I, test_input_t const &tsi)
	{
		// the reference point
		x_t x0 = tsi.get_x();

		// perform inverse mapping to obtain image of reference point on element
		double tol = 1e-6;
		unsigned max_iter = 100;
		inverse_mapping<elem_t> im(m_elem);
		if (!im.eval(x0, tol, max_iter))
			throw std::runtime_error("Could not perform inverse mapping");
		auto res = im.get_result();
		m_xi0 = res.topRows(xi_t::RowsAtCompileTime);
		m_zeta = res(xi_t::RowsAtCompileTime);

		// perform telles transform
		quad_t quad_telles = telles_transform(m_quad_orig, m_xi0, m_zeta);

		for (auto q : quad_telles)
		{
			xi_t xi = q.get_xi();
			double w = q.get_w();
			w_trial_input_t wtri(m_elem, xi);
			auto N = trial_nset_t::template eval_shape<0>(xi);
			I += semi_block_product(m_kernel(tsi, wtri) * wtri.get_jacobian() * w, N);
		}

	} // end of function integrate

private:
	elem_t const &m_elem;		/**< \brief the element reference */
	kernel_base<kernel_t> const &m_kernel;	/**< \brief the kernel reference */
	xi_t m_xi0;
	double m_zeta;
	quad_t m_quad_orig;
};

} // end of namespace NiHu

#endif /* NIHU_NEARLY_SINGULAR_COLLOCATIONAL_HPP_INCLUDED */
