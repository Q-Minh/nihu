/**
 * \file unit_kernel.hpp
 * \ingroup library
 * \brief implementation of the unit kernel \f$K(x,y) = 1\f$
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef UNIT_KERNEL_HPP_INCLUDED
#define UNIT_KERNEL_HPP_INCLUDED

#include "../bem/kernel.hpp"
#include "../bem/gaussian_quadrature.hpp"

#include "empty_input.hpp"

/** \brief the unit kernel returning K(x,y) = 1 for all inputs */
template <class Space>
class unit_kernel;

/** \brief traits of the 3D unit kernel */
template<class Space>
struct kernel_traits<unit_kernel<Space> >
{
	/** \brief test input type */
	typedef empty_input<Space> test_input_t;
	/** \brief trial input type */
	typedef empty_input<Space> trial_input_t;
	/** \brief kernel result type */
	typedef typename Space::scalar_t result_t;
	/** \brief quadrature family tag */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief shows if kernel is symmetric */
	static bool const is_symmetric = true;
	/** \brief the singularity order (r^-order) */
	static unsigned const singularity_order = 0;
	/** \brief the singular quadrature order */
	static unsigned const singular_quadrature_order = 0;
};

/** \brief the unit kernel returning K(x,y) = 1 for all inputs */
template <class Space>
class unit_kernel :
	public kernel_base<unit_kernel<Space> >
{
public:
	typedef kernel_base<unit_kernel<Space> > base_t;
	typedef typename base_t::scalar_t scalar_t;
	typedef typename base_t::result_t result_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;

	/**
	 * \brief evaluate kernel at test and trial positions
	 * \param [in] x test position
	 * \param [in] y trial position
	 * \return kernel value K(x,y)
	 */
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		return 1.0;
	}
	
	/**
	 * \brief evaluate kernel complexity at
	 * \param [in] x test position
	 * \param [in] y trial position
	 * \param [in] s linear element size
	 * \return kernel value K(x,y)
	 */
	unsigned estimate_complexity(
		test_input_t const &x,
		trial_input_t const &y,
		scalar_t const &s) const
	{
		return 0;
	}
};

#endif // UNIT_KERNEL_HPP_INCLUDED

