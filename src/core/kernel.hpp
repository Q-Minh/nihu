// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2013  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2013  Peter Rucz <rucz@hit.bme.hu>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

/**
 * \file kernel.hpp
 * \ingroup kernel
 * \brief implementation of various kernels
 */

#ifndef KERNEL_HPP_INCLUDED
#define KERNEL_HPP_INCLUDED

#include "complexity_estimator.hpp"
#include "../util/crtp_base.hpp"
#include "../util/brick.hpp"
#include "../util/collection.hpp"
#include "../util/couple.hpp"

class empty_data {};

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
	NIHU_CRTP_HELPERS

	/** \brief the traits class */
	typedef kernel_traits<Derived> traits_t;

	/** \brief type of the first (test) kernel input */
	typedef typename traits_t::test_input_t test_input_t;
	/** \brief type of the second (trial) kernel input */
	typedef typename traits_t::trial_input_t trial_input_t;
	/** \brief the kernel data type */
	typedef typename traits_t::data_t data_t;
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

	/** \brief constructor initialising the kernel data */
	kernel_base(data_t const &data = data_t()) :
		m_data(data)
	{
	}

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
		output_t output(x, y, m_data);
		// select result from output
		return output.get_result();
	}

	data_t const &get_data(void) const
	{
		return m_data;
	}

private:
	/** \brief the kernel data */
	data_t m_data;
};


/** \brief a set of kernel outputs gathered in a couple
 * \tparam outputs the output classes
 */
template <class...outputs>
class couple_output :
	public merge<outputs...>::type
{
public:
	typedef typename merge<outputs...>::type base_t;
	typedef couple<typename outputs::result_t...> result_t;

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
		return result_t(
			static_cast<typename find_in_wall<outputs, base_t>::type const &>(*this).get_result()...
		);
	}
};


// forward declaration
template <class...Kernels>
class couple_kernel;

/** \brief specialisation of ::kernel_traits for the ::couple_kernel class */
template <class...Kernels>
struct kernel_traits<couple_kernel<Kernels...> >
{
	/** \brief type of the first (test) kernel input */
	typedef typename merge<typename kernel_traits<Kernels>::test_input_t...>::type test_input_t;
	/** \brief type of the second (trial) kernel input */
	typedef typename merge<typename kernel_traits<Kernels>::trial_input_t...>::type trial_input_t;
	/** \brief the data type */
	typedef typename merger<typename Kernels::data_t...>::ret_type data_t;
	/** \brief type of the kernel output (not the result) */
	typedef couple_output<typename kernel_traits<Kernels>::output_t...> output_t;
	/** \brief type of the kernel's result */
	typedef typename output_t::result_t result_t;
	/** \brief the quadrature family the kernel is integrated with
	\todo static assert here or something more clever
	*/
	typedef typename kernel_traits<
		typename std::tuple_element<0, std::tuple<Kernels...> >::type
	>::quadrature_family_t quadrature_family_t;
	/** \brief true if K(x,y) = K(y,x) */
	static bool const is_symmetric = tmp::and_<
		std::integral_constant<bool, kernel_traits<Kernels>::is_symmetric>...
	>::value;
	/** \brief the singularity order ( r^(-order) ) */
	static unsigned const singularity_order = tmp::max_<
		std::integral_constant<unsigned, kernel_traits<Kernels>::singularity_order>...
	>::value;
	/** \brief the quadrature order used for the generation of Duffy type singular quadratures */
	static unsigned const singular_quadrature_order = tmp::max_<
		std::integral_constant<unsigned, kernel_traits<Kernels>::singular_quadrature_order>...
	>::value;
	/** \brief the kernel complexity estimator class */
	typedef typename merge_kernel_complexity_estimators<
		typename kernel_traits<Kernels>::complexity_estimator_t...
	>::type complexity_estimator_t;
};


/** \brief a kernel consisting of a set of kernels
 * \tparam Kernels the type list of kernels
 */
template <class...Kernels>
class couple_kernel :
	public kernel_base<couple_kernel<Kernels...> >
{
public:
	/** \brief the traits class */
	typedef kernel_base<couple_kernel<Kernels...> > base_t;

	/** \brief constructor from list of constant references
	 * \param [in] kernels references to kernel instances
	 */
	couple_kernel(kernel_base<Kernels> const &...kernels) :
		base_t(merge_data(kernels.get_data()...))
	{
	}
};


/** \brief factory function to create a couple kernel instance
 * \tparam Args type list of kernels
 * \param [in] args the kernel instances
 * \return a couple kernel instance
 */
template <class...Args>
couple_kernel<Args...>
	create_couple_kernel(kernel_base<Args> const &...args)
{
	return couple_kernel<Args...>(args...);
}

#endif // KERNEL_HPP_INCLUDED

