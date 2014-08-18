// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
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
 * \brief implementation of class kernel and its traits, as well as couple kernels
 */

#ifndef KERNEL_HPP_INCLUDED
#define KERNEL_HPP_INCLUDED

#include "complexity_estimator.hpp"
#include "asymptotic_types.hpp"
#include "../util/crtp_base.hpp"
#include "../util/brick.hpp"
#include "../util/collection.hpp"
#include "../util/couple.hpp"
#include "../library/interval_estimator.hpp"
#include "../library/distance_kernel_intervals.hpp"

class empty_data {};

/** \brief metafunctions returning regular and singular kernel traits */
namespace kernel_traits_ns
{
	/** \brief return the coordinate space where the kernel is defined */
	template <class Derived> struct space;
	/** \brief return the kernel's test input */
	template <class Derived> struct test_input;
	/** \brief return the kernel's trial input */
	template <class Derived> struct trial_input;
	/** \brief return the kernel's data type */
	template <class Derived> struct data;
	/** \brief return the kernel's output type */
	template <class Derived> struct output;
	/** \brief return the quadrature family the kerenel is integrated with */
	template <class Derived> struct quadrature_family;
	/** \brief return the far field asymptotic behaviour of the kernel */
	template <class Derived> struct far_field_behaviour;

	/** \brief return the kernel's result dimensionality
	 * \todo this should be automatically deduced from the kernel result
	 */
	template <class Derived> struct result_dimension;
	/** \brief return whether the kernel is symmetric or not */
	template <class Derived> struct is_symmetric;
	/** \brief return whether the kernel is singular or not */
	template <class Derived> struct is_singular;

	/** \brief return the kernel's singularity type */
	template <class Derived> struct singularity_type;
	/** \brief return the quadrature order the singular kernel needs to be integrated with */
	template <class Derived> struct singular_quadrature_order;
	/** \brief return the kernel's singular core type */
	template <class Derived> struct singular_core;
}

/**
* \brief traits class of a kernel
* \tparam Derived the CRTP derived kernel
*/
template <class Derived>
struct kernel_traits
{
	typedef typename kernel_traits_ns::space<Derived>::type space_t;
	/** \brief kernel test input type */
	typedef typename kernel_traits_ns::test_input<Derived>::type test_input_t;
	/** \brief kernel trial input type */
	typedef typename kernel_traits_ns::trial_input<Derived>::type trial_input_t;
	/** \brief the data type */
	typedef typename kernel_traits_ns::data<Derived>::type data_t;
	/** \brief the kernel output type */
	typedef typename kernel_traits_ns::output<Derived>::type output_t;
	/** \brief the far field asymptotic behaviour of the kernel */
	typedef typename kernel_traits_ns::far_field_behaviour<Derived>::type far_field_behaviour_t;
	/** \brief the far field asymptotic behaviour of the kernel */
	typedef typename kernel_traits_ns::quadrature_family<Derived>::type quadrature_family_t;



	enum {
		/** \brief the kernel result's dimension
		 * \todo should be computed from the kernel result and not defined separately in traits
		 */
		result_dimension = kernel_traits_ns::result_dimension<Derived>::value,
		/** \brief indicates if the kernel is symmetric */
		is_symmetric = kernel_traits_ns::is_symmetric<Derived>::value,
		/** \brief indicates if the kernel is singular */
		is_singular = kernel_traits_ns::is_singular<Derived>::value
	};
};


template <class Derived>
struct singular_kernel_traits
{
	typedef typename kernel_traits_ns::singularity_type<Derived>::type singularity_type_t;

	typedef typename kernel_traits_ns::singular_core<Derived>::type singular_core_t;

	enum {
		singular_quadrature_order = kernel_traits_ns::singular_quadrature_order<Derived>::value
	};
};



template <class Derived>
struct kernel_compl_estimator
{
	typedef interval_estimator<
		typename distance_kernel_interval<
			typename kernel_traits<Derived>::far_field_behaviour_t
		>::type
	> type;
};


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
	/** \brief compile time check if the two kernel inputs are compatible */
	static_assert(std::is_same<typename test_input_t::space_t, typename trial_input_t::space_t>::value,
		"The test and trial kernel inputs must define the same coordinate space");
	/** \brief the kernel data type */
	typedef typename traits_t::data_t data_t;
	/** \brief type of the kernel output (not the result) */
	typedef typename traits_t::output_t output_t;
	/** \brief kernel result type */
	typedef typename output_t::result_t result_t;
	/** \brief the kernel result's dimensionality */
	enum { result_dimension = traits_t::result_dimension };

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

	/** \brief true if the kernel is singular */
	static bool const is_sungular = traits_t::is_singular;

	/** \brief the asymptotic (far field) behaviour of the kernel */
	typedef typename traits_t::far_field_behaviour_t far_field_behaviour_t;

	/** \brief the kernel complexity estimator class */
	typedef typename kernel_compl_estimator<Derived>::type estimator_t;

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

	/** \brief return kernel data
	 * \return kernel data
	 */
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
private:
	/** \brief this metafunction is an msvc workaround for variadic templates */
	template <class T>
	struct result_mf
	{
		typedef typename T::result_t type;
	};

public:
	/** \brief the base type */
	typedef typename merge<outputs...>::type base_t;
	/** \brief the merged output type */
	typedef couple<typename result_mf<outputs>::type...> result_t;

	/** \brief constructor
	 * \tparam test_input_t the test input type
	 * \tparam trial_input_t the trial input type
	 * \tparam kernel_t the kernel type
	 * \param [in] test_input the test input
	 * \param [in] trial_input the trial input
	 * qparam [in] kernel the kernel instance
	 */
	template <class test_input_t, class trial_input_t, class kernel_t>
	couple_output(
		test_input_t const &test_input,
		trial_input_t const &trial_input,
		kernel_t const &kernel) :
		base_t(test_input, trial_input, kernel)
	{
	}

	/** \brief return result
	 * \return couple result
	 */
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
public:
	/** \brief type of the first (test) kernel input */
	typedef typename merge<typename kernel_traits<Kernels>::test_input_t...>::type test_input_t;
	/** \brief type of the second (trial) kernel input */
	typedef typename merge<typename kernel_traits<Kernels>::trial_input_t...>::type trial_input_t;
	/** \brief the data type */
	typedef typename merger<typename kernel_traits<Kernels>::data_t...>::ret_type data_t;
	/** \brief type of the kernel output (not the result) */
	typedef couple_output<typename kernel_traits<Kernels>::output_t...> output_t;
	/** \brief type of the kernel's result */
	typedef typename output_t::result_t result_t;
	/** \brief the kernel result's dimensionality
	 * \todo update this, couple does not work with vector kernels with this definition
	 */
	enum { result_dimension = 1 };
	/** \brief the quadrature family the kernel is integrated with
	 * \todo static assert here or something more clever
	 */
	typedef typename kernel_traits<
		typename std::tuple_element<0, std::tuple<Kernels...> >::type
	>::quadrature_family_t quadrature_family_t;
	/** \brief true if K(x,y) = K(y,x) */
	static bool const is_symmetric = tmp::and_<
		std::integral_constant<bool, kernel_traits<Kernels>::is_symmetric>...
	>::value;
	/** \brief true if any of the kernels is singular */
	static bool const is_singular = tmp::or_<
		std::integral_constant<bool, kernel_traits<Kernels>::is_singular>...
	>::value;
	/** \brief the combined far field behaviour order */
	typedef typename tmp::max_<
		typename kernel_traits<Kernels>::far_field_behaviour_t...
	>::type far_field_behaviour_t;
};

template <class...Kernels>
struct kernel_compl_estimator<couple_kernel<Kernels...> >
{
	typedef typename merge_kernel_complexity_estimators<
		typename kernel_compl_estimator<Kernels>::type...
	>::type type;
};

/** \brief specialisation of ::singular_kernel_traits for the ::couple_kernel class */
template <class...Kernels>
struct singular_kernel_traits<couple_kernel<Kernels...> >
{
private:
	template <class K>
	struct sing_quad_order_constant : tmp::integer<unsigned, singular_kernel_traits<K>::singular_quadrature_order> {};

public:
	/** \brief the combined singularity order */
	typedef typename tmp::max_<
		typename singular_kernel_traits<Kernels>::singularity_type_t...
	>::type singularity_type_t;
	/** \brief the quadrature order used for the generation of blind singular quadratures */
	static unsigned const singular_quadrature_order = tmp::max_<
		typename sing_quad_order_constant<Kernels>::type...
	>::type::value;
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
couple_kernel<Args...> create_couple_kernel(kernel_base<Args> const &...args)
{
	return couple_kernel<Args...>(args...);
}

#endif // KERNEL_HPP_INCLUDED

