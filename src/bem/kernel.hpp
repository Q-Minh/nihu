/**
 * \file kernel.hpp
 * \brief implementation of various kernels
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */
#ifndef KERNEL_HPP_INCLUDED
#define KERNEL_HPP_INCLUDED

#include "kernel_input.hpp"
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

	/** \compile time check if the two kernel inputs are compatible */
	static_assert(std::is_same<typename test_input_t::space_t, typename trial_input_t::space_t>::value,
		"The test and trial kernel inputs must define the same coordinate space");

	/** \brief type of the kernel's domain space */
	typedef typename test_input_t::space_t space_t;
	/** \brief type of a location vector in the kernel's domain */
	typedef typename space_t::location_t x_t;
	/** \brief type of the scalar coordinate in the kernel's domain */
	typedef typename space_t::scalar_t scalar_t;

	/** \brief type of the kernel's result */
	typedef typename traits_t::result_t result_t;
	
	/** \brief the quadrature family the kernel is integrated with */
	typedef typename traits_t::quadrature_family_t quadrature_family_t;
	
	/** \brief true if K(x,y) = K(y,x) */
	static bool const is_symmetric = traits_t::is_symmetric;
	/** \brief the singularity order ( r^(-order) ) */
	static unsigned const singularity_order = traits_t::singularity_order;
	/** \brief the quadrature order used for the generation of Duffy type singular quadratures */
	static unsigned const singular_quadrature_order = traits_t::singular_quadrature_order;

private:
	/** \brief CRTP helper function */
	Derived const &derived(void) const
	{
		return static_cast<Derived const &>(*this);
	}

	/** \brief CRTP helper function */
	Derived &derived(void)
	{
		return static_cast<Derived &>(*this);
	}

public:
	/** \brief the kernel bound at the test kernel input */
	class kernel_bind
	{
	public:
		/** \brief constructor
		* \param [in] kernel the kernel that is bound
		* \param [in] x the test input that binds the kernel
		*/
		kernel_bind(Derived &kernel, test_input_t const &x)
			: m_kernel(kernel), m_test_input(x)
		{
		}

		/** \brief evaluate bound kernel
		* \param [in] y the trial input
		* \return the kernel result value
		*/
		result_t eval(trial_input_t const &y) const
		{
			return m_kernel.eval(m_test_input, y);
		}

	private:
		/** \brief the kernel that is bound */
		Derived &m_kernel;
		/** \brief the test input used to bind the kernel */
		test_input_t const &m_test_input; 
	};

	/** \brief constructor */
	kernel_base() :
		m_num_evaluations(0)
	{
	}

	/** \brief bind the kernel at its test input
	* \param [in] x the test input that binds the kernel
	* \return the bound kernel function
	*/
	kernel_bind bind(test_input_t const &x)
	{
		return kernel_bind(derived(), x);
	}

	/**
	 * \brief return number of kernel evaluations
	 * \return number of kernel evaluations
	 */
	long long unsigned get_num_evaluations(void) const
	{
		return m_num_evaluations;
	}

	/**
	 * \brief evaluate kernel at a given source and receiver position
	 * \param [in] x test position
	 * \param [in] x trial position
	 * \return kernel value K(x,y)
	 */
	result_t eval(test_input_t const &x, trial_input_t const &y)
	{
		m_num_evaluations++;
		return derived()(x, y);
	}

	/**
	 * \brief determine kernel's polynomial complexity
	 * \param [in] x test position
	 * \param [in] y trial position
	 * \param [in] reference_size linear estimated size of the trial element
	 * \return polynomial degree needed for accurate integration
	 */
	unsigned estimate_complexity(test_input_t const &x, trial_input_t const &y, scalar_t const &reference_size) const
	{
		return derived().estimate_complexity(x, y, reference_size);
	}

protected:
	/** \brief number of kernel evaluations */
	long long unsigned m_num_evaluations;
};


#endif // KERNEL_HPP_INCLUDED

