/**
 * \file unit_kernel.hpp
 * \ingroup library
 * \brief implementation of the unit kernel \f$K(x,y) = 1\f$
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef UNIT_KERNEL_HPP_INCLUDED
#define UNIT_KERNEL_HPP_INCLUDED

#include "../bem/kernel.hpp"
#include "../util/brick.hpp"
#include "../bem/gaussian_quadrature.hpp"

/** \brief the unit kernel returning K(x,y) = 1 for all inputs */
template <class Scalar>
class unit_kernel;

/** \brief traits of the unit kernel */
template<class Scalar>
struct kernel_traits<unit_kernel<Scalar> >
{
	/** \brief test input type */
	typedef empty_wall test_input_t;
	/** \brief trial input type */
	typedef empty_wall trial_input_t;
	/** \brief kernel result type */
	typedef Scalar result_t;
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
template <class Scalar>
class unit_kernel :
	public kernel_base<unit_kernel<Scalar> >
{
public:
	/** \brief the crtp base type */
	typedef kernel_base<unit_kernel<Scalar> > base_t;
	/** \brief the scalar type */
	typedef typename base_t::scalar_t scalar_t;
	/** \brief the result type */
	typedef typename base_t::result_t result_t;
	/** \brief the test input type */
	typedef typename base_t::test_input_t test_input_t;
	/** \brief the trial input type */
	typedef typename base_t::trial_input_t trial_input_t;

	/**
	 * \brief evaluate kernel at test and trial positions
	 * \param [in] x test position
	 * \param [in] y trial position
	 * \return kernel value K(x,y)
	 */
	constexpr result_t operator()(
		test_input_t const &, trial_input_t const &) const
	{
		return result_t(1.0);
	}
	
	/**
	 * \brief evaluate kernel complexity at
	 * \param [in] x test position
	 * \param [in] y trial position
	 * \param [in] s linear element size
	 * \return kernel value K(x,y)
	 */
	constexpr unsigned estimate_complexity(
		test_input_t const &, trial_input_t const &, scalar_t const &) const
	{
		return 0;
	}
};

#endif // UNIT_KERNEL_HPP_INCLUDED

