/**
 * \file kernel.hpp
 * \ingroup intop
 * \brief implementation of various kernels
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */
 
#ifndef KERNEL_HPP_INCLUDED
#define KERNEL_HPP_INCLUDED

#include "../util/crtp_base.hpp"
#include "kernel_input.hpp"
#include "../util/brick.hpp"
#include "couple.hpp"

/**
* \brief traits class of a kernel
* \tparam Derived the CRTP derived kernel
*/
template <class Derived>
struct kernel_traits;

/**
 * \brief CRTP base class of all BEM kernels
 * \tparam Derived the CRTP derived class
 */
template <class Derived>
class kernel_base
{
public:
	/** \brief the traits class */
	typedef kernel_traits<Derived> traits_t;
	
	/** \brief type of the first (test) kernel input */
	typedef typename traits_t::test_input_t test_input_t;
	/** \brief type of the second (trial) kernel input */
	typedef typename traits_t::trial_input_t trial_input_t;
	/** \brief type of the kernel output (not the result) */
	typedef typename traits_t::output_t output_t;
	/** \brief type of the kernel's result */
	typedef typename traits_t::result_t result_t;

	/** \brief compile time check if the two kernel inputs are compatible */
	static_assert(std::is_same<typename test_input_t::space_t, typename trial_input_t::space_t>::value,
		"The test and trial kernel inputs must define the same coordinate space");

	/** \brief type of the kernel's domain space */
	typedef typename test_input_t::space_t space_t;
	/** \brief type of a location vector in the kernel's domain */
	typedef typename space_t::location_t x_t;
	/** \brief type of the scalar coordinate in the kernel's domain */
	typedef typename space_t::scalar_t scalar_t;
	
	/** \brief the quadrature family the kernel is integrated with */
	typedef typename traits_t::quadrature_family_t quadrature_family_t;
	
	/** \brief true if K(x,y) = K(y,x) */
	static bool const is_symmetric = traits_t::is_symmetric;
	/** \brief the singularity order ( r^(-order) ) */
	static unsigned const singularity_order = traits_t::singularity_order;
	/** \brief the quadrature order used for the generation of Duffy type singular quadratures */
	static unsigned const singular_quadrature_order = traits_t::singular_quadrature_order;

private:
	NIHU_CRTP_HELPERS

public:
	/** \brief the kernel bound at the test kernel input */
	class kernel_bind
	{
	public:
		/** \brief constructor
		* \param [in] kernel the kernel that is bound
		* \param [in] x the test input that binds the kernel
		*/
		kernel_bind(kernel_base<Derived> const &kernel, test_input_t const &x)
			: m_kernel(kernel.derived()), m_test_input(x)
		{
		}

		/** \brief evaluate bound kernel
		* \param [in] y the trial input
		* \return the kernel result value
		*/
		result_t operator()(trial_input_t const &y) const
		{
			return m_kernel(m_test_input, y);
		}

	private:
		/** \brief the kernel that is bound */
		Derived const &m_kernel;
		/** \brief the test input used to bind the kernel */
		test_input_t const &m_test_input; 
	};

	/** \brief bind the kernel at its test input
	* \param [in] x the test input that binds the kernel
	* \return the bound kernel function
	*/
	kernel_bind bind(test_input_t const &x) const
	{
		return kernel_bind(*this, x);
	}

	/**
	 * \brief evaluate kernel at a given source and receiver position
	 * \param [in] x test position
	 * \param [in] y trial position
	 * \return kernel value K(x,y)
	 */
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		// instantiate output
		output_t output(x, y, *this);
		// select result from output
		return output.get_result();
	}

	/**
	 * \brief determine kernel's polynomial complexity
	 * \param [in] x test position
	 * \param [in] y trial position
	 * \param [in] reference_size linear estimated size of the trial element
	 * \return polynomial degree needed for accurate integration
	 */
	unsigned estimate_complexity_interface(
		test_input_t const &x,
		trial_input_t const &y,
		scalar_t const &reference_size) const
	{
		return derived().estimate_complexity(x, y, reference_size);
	}
};


template <class out1, class out2>
class couple_output :
	merge<out1, out2>::type
{
public:
	typedef typename merge<out1, out2>::type base_t;
	typedef couple<typename out1::result_t, typename out2::result_t> result_t;

	template <class test_input_t, class trial_input_t, class kernel_t>
	couple_output(
		test_input_t const &test_input,
		trial_input_t const &trial_input,
		kernel_t const &kernel) :
		base_t(test_input, trial_input, kernel)
	{
	}

	result_t get_result(void) const
	{
		return create_couple(
			static_cast<typename find_in_wall<out1, base_t>::type const &>(*this).get_result(),
			static_cast<typename find_in_wall<out2, base_t>::type const &>(*this).get_result()
		);
	}
};



template <class Kernel1, class Kerel2>
class couple_kernel;

template <class Kernel1, class Kernel2>
struct kernel_traits<couple_kernel<Kernel1, Kernel2> >
{
	/** \brief type of the first (test) kernel input */
	typedef typename merge<
		typename kernel_traits<Kernel1>::test_input_t,
		typename kernel_traits<Kernel2>::test_input_t
	>::type test_input_t;
	/** \brief type of the second (trial) kernel input */
	typedef typename merge<
		typename kernel_traits<Kernel1>::trial_input_t,
		typename kernel_traits<Kernel2>::trial_input_t
	>::type trial_input_t;
	/** \brief type of the kernel output (not the result) */
	typedef couple_output<
		typename kernel_traits<Kernel1>::output_t,
		typename kernel_traits<Kernel2>::output_t
	> output_t;
	/** \brief type of the kernel's result */
	typedef typename output_t::result_t result_t;
	/** \brief the quadrature family the kernel is integrated with
	\todo static assert here or something more clever
	*/
	typedef typename kernel_traits<Kernel1>::quadrature_family_t quadrature_family_t;
	/** \brief true if K(x,y) = K(y,x) */
	static bool const is_symmetric =
		kernel_traits<Kernel1>::is_symmetric &&
		kernel_traits<Kernel2>::is_symmetric ;
	/** \brief the singularity order ( r^(-order) ) */
	static unsigned const singularity_order = std::max(
		kernel_traits<Kernel1>::singularity_order,
		kernel_traits<Kernel1>::singularity_order);
	/** \brief the quadrature order used for the generation of Duffy type singular quadratures */
	static unsigned const singular_quadrature_order = std::max(
		kernel_traits<Kernel1>::singular_quadrature_order,
		kernel_traits<Kernel1>::singular_quadrature_order);
};


template <class Kernel1, class Kernel2>
class couple_kernel :
	public kernel_base<couple_kernel<Kernel1, Kernel2> >
{
public:
	couple_kernel(
		kernel_base<Kernel1> const &k1,
		kernel_base<Kernel2> const &k2)
	{
	}
};


template <class Kernel1, class Kernel2>
couple_kernel<Kernel1, Kernel2>
	create_couple_kernel(kernel_base<Kernel1> const &k1, kernel_base<Kernel2> const &k2)
{
	return couple_kernel<Kernel1, Kernel2>(k1, k2);
}

#endif // KERNEL_HPP_INCLUDED

