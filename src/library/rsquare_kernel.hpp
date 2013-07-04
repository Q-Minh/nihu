/**
 * \file rsquare_kernel.hpp
 * \ingroup library
 * \brief implementation of kernel \f$1/\sqrt{r}\f$
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef RSQUARE_KERNEL_HPP_INCLUDED
#define RSQUARE_KERNEL_HPP_INCLUDED

#include "../bem/kernel.hpp"
#include "../bem/gaussian_quadrature.hpp"

#include "location_and_normal.hpp"
#include "reciprocal_distance_kernel.hpp"

// forward declaration
class rsquare_kernel;

/** \brief traits of the poisson G kernel */
template<>
struct kernel_traits<rsquare_kernel>
{
	/** \brief kernel test input type */
	typedef location<space_2d> test_input_t;
	/** \brief kernel trial input type */
	typedef location<space_2d> trial_input_t;
	/** \brief kernel result type */
	typedef space_2d::scalar_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 1;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
};

/** \brief 3D Poisson kernel \f$1/4\pi r\f$ */
class rsquare_kernel :
	public kernel_base<rsquare_kernel>,
	public reciprocal_distance_kernel<rsquare_kernel>
{
public:
	/** \brief the crtp base's type */
	typedef kernel_base<rsquare_kernel> base_t;
	/** \brief type of the test input */
	typedef typename base_t::test_input_t test_input_t;
	/** \brief type of the trial input */
	typedef typename base_t::trial_input_t trial_input_t;
	/** \brief type of the scalar of the inputs */
	typedef typename base_t::scalar_t scalar_t;

	/** \brief evaluate kernel between test and trial positions
	* \param [in] x the test input
	* \param [in] y the trial input
	* \return the kernel value K(x,y)
	*/
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		scalar_t r = (y.get_x() - x.get_x()).norm();
		return 1.0 / std::sqrt(r);
	}
};



#endif // RSQUARE_KERNEL_HPP_INCLUDED


