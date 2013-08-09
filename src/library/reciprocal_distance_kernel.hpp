/**
 * \file reciprocal_distance_kernel.hpp
 * \ingroup library
 * \brief complexity rules for kernels with a \f$ 1/r^n\f$ behaviour
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef RECIPROCAL_DISTANCE_KERNEL
#define RECIPROCAL_DISTANCE_KERNEL

/** \brief class to estimate polynomial complexity of \f$1/r^n\f$ type kernels */
template <class Derived>
class reciprocal_distance_kernel
{
private:
	/** \brief the kernel traits */
	typedef kernel_traits<Derived> traits_t;
	/** \brief the test input */
	typedef typename traits_t::test_input_t test_input_t;
	/** \brief the trial input */
	typedef typename traits_t::trial_input_t trial_input_t;
	/** \brief the kernel's scalar type */
	typedef typename test_input_t::scalar_t scalar_t;

public:
	/** \brief estimate polynomial complexity
	* \param [in] x test input
	* \param [in] y trial input
	* \param [in] s relative distance
	* \return polynomial complexity
	*/
	unsigned estimate_complexity(
		test_input_t const &x,
		trial_input_t const &y,
		scalar_t const &s) const
	{
		return estimate_complexity_impl(
			std::integral_constant<unsigned, traits_t::singularity_order>(),
			(x.get_x() - y.get_x()).norm() / s);
	}

private:
	unsigned estimate_complexity_impl(
		std::integral_constant<unsigned, 1>,
		scalar_t const &rel_dist) const
	{
		/** \todo find a better approach */
		return rel_dist < 2.0 ? 7 : ( rel_dist < 5.0 ? 3 : 0);
	}

	unsigned estimate_complexity_impl(
		std::integral_constant<unsigned, 2>,
		scalar_t const &rel_dist) const
	{
		/** \todo find a better approach */
		return rel_dist < 2.0 ? 5 : ( rel_dist < 5.0 ? 3 : 0);
	}

	unsigned estimate_complexity_impl(
		std::integral_constant<unsigned, 3>,
		scalar_t const &rel_dist) const
	{
		/** \todo This is a blind copy, find a better approach! */
		return rel_dist < 2.0 ? 5 : ( rel_dist < 5.0 ? 3 : 0);
	}
};

#endif // RECIPROCAL_DISTANCE_KERNEL

