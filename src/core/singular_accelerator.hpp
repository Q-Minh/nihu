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

/// \file singular_accelerator.hpp
/// \ingroup quadrature
/// \brief implementation of class NiHu::singular_accelerator

#ifndef SINGULAR_ACCELERATOR_HPP_INCLUDED
#define SINGULAR_ACCELERATOR_HPP_INCLUDED

#include "asymptotic_types.hpp"
#include "kernel.hpp"	// for the galerkin case
#include "singular_galerkin_quadrature.hpp"	// for the galerkin case
#include "blind_transform_selector.hpp"	// for the collocation case
#include "blind_singular_quadrature.hpp"	// for the collocation case
#include "field.hpp"
#include "../tmp/control.hpp"
#include "../util/is_specialisation.hpp"

#include "singular_integral_shortcut.hpp"

#include "../util/dual_range.hpp"
#include "formalism.hpp"

#include <stdexcept>

namespace NiHu
{

/**
* \brief a dual iterator to point to a test and a trial quadrature element
* \tparam test_iterator_t the iterator type of the test quadrature
* \tparam trial_iterator_t the iterator type of the trial quadrature
*/
template <class test_iterator_t, class trial_iterator_t>
class singular_quadrature_iterator :
	public dual_iterator<iteration::diagonal, test_iterator_t, trial_iterator_t>
{
public:
	/** \brief the base type */
	typedef dual_iterator<iteration::diagonal, test_iterator_t, trial_iterator_t> base_t;
	/** \brief the value type of both quadratures */
	typedef typename test_iterator_t::value_type value_type;

	/** \brief constructor initialising members
	* \param [in] test the test quadrature
	* \param [in] trial the trial quadrature
	*/
	singular_quadrature_iterator(test_iterator_t test, trial_iterator_t trial)
		: base_t(test, trial)
	{
	}

	/** \brief return the test quadrature element
	* \return test quadrature element
	*/
	value_type const &get_test_quadrature_elem(void) const
	{
		return *base_t::get_first();
	}

	/** \brief return the trial quadrature element
	* \return trial quadrature element
	*/
	value_type const &get_trial_quadrature_elem(void) const
	{
		return *base_t::get_second();
	}
};


/** \brief an invalid singular iterator assigned to nonexisting integrals */
class invalid_singular_iterator {};

/** \brief an invalid singular accelerator assigned to nonexisting integrals */
class invalid_singular_accelerator
{
public:
	/** \brief return the begin iterator
	 * \tparam match the element match type
	 * \return returns and invalud iterator
	 */
	template <class match>
	invalid_singular_iterator begin(match const &) const
	{
		return invalid_singular_iterator();
	}

	/** \brief return the end iterator
	 * \tparam match the element match type
	 * \return returns and invalud iterator
	 */
	template <class match>
	invalid_singular_iterator end(match const &) const
	{
		return invalid_singular_iterator();
	}
};

/**
* \brief an accelerator class that stores singular quadratures for different singularity types
* \tparam Kernel the kernel that is integrated
* \tparam TestField the test field over which integration is performed
* \tparam TrialField the trial field over which integration is performed
*/
template <
	class Kernel,
	class TestField,
	class TrialField,
	class = typename get_formalism<TestField, TrialField>::type
>
class singular_accelerator;


/**
* \brief specialisation of NiHu::singular_accelerator for the general formalism
* \tparam Kernel the kernel that is integrated
* \tparam TestField the test field over which integration is performed
* \tparam TrialField the trial field over which integration is performed
*/
template <
	class Kernel,
	class TestField,
	class TrialField
>
class singular_accelerator<Kernel, TestField, TrialField, formalism::general>
{
	// CRTP check
	static_assert(std::is_base_of<kernel_base<Kernel>, Kernel>::value,
		"The kernel field must be derived from kernel_base<Kernel>");
	static_assert(std::is_base_of<field_base<TestField>, TestField>::value,
		"The test field must be derived from field_base<TestField>");
	static_assert(std::is_base_of<field_base<TrialField>, TrialField>::value,
		"The trial field must be derived from field_base<TrialField>");
public:
	/** \brief template argument as nested type */
	typedef Kernel kernel_t;
	/** \brief template argument as nested type */
	typedef TestField test_field_t;
	/** \brief template argument as nested type */
	typedef TrialField trial_field_t;
	
	/** \brief the test eleme type */
	typedef typename test_field_t::elem_t test_elem_t;
	
	/** \brief indicates if the elements are surface elements or not */
	static bool const is_surface = element_traits::is_surface_element<test_elem_t>::type::value;

	/** \brief the L-set of the test domain */
	typedef typename test_field_t::elem_t::lset_t test_lset_t;
	/** \brief the L-set of the trial domain */
	typedef typename trial_field_t::elem_t::lset_t trial_lset_t;

	/** \brief test domain */
	typedef typename test_field_t::elem_t::domain_t test_domain_t;
	/** \brief trial domain */
	typedef typename trial_field_t::elem_t::domain_t trial_domain_t;
	
	/** \brief the domain dimension */
	static unsigned const domain_dimension = test_domain_t::dimension;
	// report compiler error if test and trial domain dimensions do not match
	static_assert(domain_dimension == trial_domain_t::dimension,
		"The test and trial domain dimensions must match");

	/** \brief number of test domain corners */
	static unsigned const n_test_corners = test_domain_t::num_corners;
	/** \brief number of trial domain corners */
	static unsigned const n_trial_corners = trial_domain_t::num_corners;

	/** \brief the scalar type of the quadrature */
	typedef typename test_domain_t::scalar_t scalar_t;
	// report compiler error if the trial and test scalar types are different
	static_assert(std::is_same<scalar_t, typename trial_domain_t::scalar_t>::value,
		"The test and trial domain scalar types must match");

	/** \brief the quadrature family obtained from the kernel */
	typedef typename kernel_traits<kernel_t>::quadrature_family_t quadrature_family_t;

	/** \brief the singular galerkin quadrature generator class */
	typedef singular_galerkin_quadrature<quadrature_family_t, test_domain_t, trial_domain_t> quad_factory_t;

	/** \brief the test quadrature type */
	typedef typename quadrature_type<quadrature_family_t, test_domain_t>::type test_quadrature_t;
	/** \brief the trial quadrature type */
	typedef typename quadrature_type<quadrature_family_t, trial_domain_t>::type trial_quadrature_t;

	/** \brief the quadrature element type (it should be the same for test and trial) */
	typedef typename test_quadrature_t::quadrature_elem_t quadrature_elem_t;
	// report compiler error if the trial and test quadrature elements are different
	static_assert(std::is_same<quadrature_elem_t, typename trial_quadrature_t::quadrature_elem_t>::value,
		"The trial and test quadrature elements must be of the same type");

	/** \brief the dual iterator type of the singular quadrature */
	typedef singular_quadrature_iterator<
		typename test_quadrature_t::const_iterator,
		typename trial_quadrature_t::const_iterator> iterator;

	/** \brief indicates whether face match is possible */
	static const bool face_match_possible = std::is_same<test_field_t, trial_field_t>::value;
	/** \brief the singular quadrature order required by the kernel */
	static unsigned const singular_quadrature_order = singular_kernel_traits<kernel_t>::singular_quadrature_order;

	/**
	 * \brief return begin iterator of the singular quadrature
	 * \return begin iterator of the singular quadrature
	 */
	iterator begin(element_match const &elem_match) const
	{
		switch (elem_match.get_match_dimension())
		{
		case domain_dimension:
			return iterator(
				m_test_quadratures[domain_dimension][0].begin(),
				m_trial_quadratures[domain_dimension][0].begin());
			break;
		case -1:
			throw std::logic_error("Cannot return singular quadrature for NO_MATCH type");
			break;
		default:
			unsigned test_domain_corner =
				test_lset_t::node_to_domain_corner(elem_match.get_overlap().get_ind1());
			unsigned trial_domain_corner =
				trial_lset_t::node_to_domain_corner(elem_match.get_overlap().get_ind2());
			return iterator(
				m_test_quadratures[elem_match.get_match_dimension()][test_domain_corner].begin(),
				m_trial_quadratures[elem_match.get_match_dimension()][trial_domain_corner].begin());
			break;
		}
	}

	/**
	* \brief return end iterator of the singular quadrature
	* \return end iterator of the singular quadrature
	*/
	iterator end(element_match const &elem_match) const
	{
		switch (elem_match.get_match_dimension())
		{
		case domain_dimension:
			return iterator(
				m_test_quadratures[domain_dimension][0].end(),
				m_trial_quadratures[domain_dimension][0].end());
			break;
		case -1:
			throw std::logic_error("Cannot return singular quadrature for NO_MATCH type");
			break;
		default:
			unsigned test_domain_corner =
				test_lset_t::node_to_domain_corner(elem_match.get_overlap().get_ind1());
			unsigned trial_domain_corner =
				trial_lset_t::node_to_domain_corner(elem_match.get_overlap().get_ind2());
			return iterator(
				m_test_quadratures[elem_match.get_match_dimension()][test_domain_corner].end(),
				m_trial_quadratures[elem_match.get_match_dimension()][trial_domain_corner].end());
			break;
		}
	}

private:
	template <class Match>
	void generate(Match, std::false_type)
	{
		int d = Match::value;
		if (d == domain_dimension) 	// face (only one)
		{
			m_test_quadratures[d].push_back(test_quadrature_t());
			m_trial_quadratures[d].push_back(trial_quadrature_t());
			quad_factory_t::template generate<Match>(
				m_test_quadratures[d][0],
				m_trial_quadratures[d][0],
				singular_quadrature_order);
		}
		else // corner or edge, lots of
		{
			test_quadrature_t test_q;
			trial_quadrature_t trial_q;
			quad_factory_t::template generate<Match>(
				test_q, trial_q, singular_quadrature_order);
				
			for	(unsigned i = 0; i < n_test_corners; ++i)
			{
				Eigen::Matrix<scalar_t, n_test_corners, domain_dimension> test_corners;
				for (unsigned k = 0; k < n_test_corners; ++k)
					test_corners.row(k) = test_domain_t::get_corner((i+k) % n_test_corners);
				m_test_quadratures[d].push_back(
					test_q.template transform<isoparam_shape_set<test_domain_t>, !is_surface>(test_corners)
					);
			}
			
			for (unsigned i = 0; i < n_trial_corners; ++i)
			{
				Eigen::Matrix<scalar_t, n_trial_corners, domain_dimension> trial_corners;
				for (unsigned k = 0; k < n_trial_corners; ++k)
				{
					unsigned idx;
					if (domain_dimension == 2 && d == 1) // quad or tria edge match
						idx = (i+1+n_trial_corners-k) % n_trial_corners;
					else
						idx = (i+k) % n_trial_corners;
					trial_corners.row(k) = trial_domain_t::get_corner(idx);
				}
				m_trial_quadratures[d].push_back(
					trial_q.template transform<isoparam_shape_set<trial_domain_t>, !is_surface>(trial_corners)
					);
			}
		}
	}
	
	template <class Match>
	void generate(Match, std::true_type)
	{
	}
	
	
private:
	template <class Match>
	struct generate_wrapper { struct type
	{
		typedef typename is_specialisation<
			singular_integral_shortcut<Kernel, TestField, TrialField, Match>
		>::type spec;
		
		void operator()(singular_accelerator &obj)
		{
			obj.generate(Match(), spec());
		}
	}; };
	
public:
	/** \brief constructor allocating and generating the quadratures */
	singular_accelerator(void)
	{
		tmp::call_each<
			typename match_type_vector<TestField, TrialField>::type,
			generate_wrapper<tmp::_1>,
			singular_accelerator &
		>(*this);
	}

private:
	std::vector<test_quadrature_t> m_test_quadratures[domain_dimension + 1];
	std::vector<trial_quadrature_t> m_trial_quadratures[domain_dimension + 1];
};


/**
 * \brief specialisation the singular accelerator for the collocational case
 * \tparam Kernel the kernel to integrate
 * \tparam TestElem the test element type
 * \tparam TrialField the trial field type
 */
template <class Kernel, class TestField, class TrialField>
class singular_accelerator<Kernel, TestField, TrialField, formalism::collocational>
{
	// CRTP check
	static_assert(std::is_base_of<kernel_base<Kernel>, Kernel>::value,
		"The kernel must be derived from kernel_base<Kernel>");
	static_assert(std::is_base_of<field_base<TrialField>, TrialField>::value,
		"The trial field must be derived from field_base<TrialField>");
public:
	/** \brief template argument as nested type */
	typedef Kernel kernel_t;
	/** \brief the test field type */
	typedef TestField test_field_t;
	/** \brief template argument as nested type */
	typedef TrialField trial_field_t;

	/** \brief the test elem type */
	typedef typename test_field_t::elem_t test_elem_t;
	/** \brief trial elem type */
	typedef typename trial_field_t::elem_t trial_elem_t;

	/** \brief trial domain type */
	typedef typename trial_elem_t::domain_t trial_domain_t;

	/** \brief trial L-set type */
	typedef typename trial_elem_t::lset_t trial_lset_t;
	/** \brief test N-set type */
	typedef typename test_field_t::nset_t test_nset_t;

	/** \brief quadrature family */
	typedef typename kernel_traits<kernel_t>::quadrature_family_t quadrature_family_t;
	/** \brief the trial quadrature type */
	typedef typename quadrature_type<
		quadrature_family_t, trial_domain_t
	>::type trial_quadrature_t;
	/** \brief quadrature element type (it should be the same for test and trial) */
	typedef typename trial_quadrature_t::quadrature_elem_t quadrature_elem_t;
	/** \brief the singular quadrature order required by the kernel */
	static unsigned const singular_quadrature_order = singular_kernel_traits<kernel_t>::singular_quadrature_order;

	/** \brief the blind transformation tag that governs the singular quadrature transformation method */
	typedef typename blind_transform_selector<
		typename singular_kernel_traits<kernel_t>::singularity_type_t,
		trial_domain_t
	>::type blind_singular_transform_tag_t;

	/** \brief the blind quadrature type */
	typedef typename blind_singular_quadrature<
		blind_singular_transform_tag_t,
		quadrature_family_t,
		trial_lset_t
	>::type trial_blind_t;

	/** \brief iterator type of the singular quadrature */
	typedef typename trial_quadrature_t::const_iterator iterator;


	/** \brief indicates whether FACE_MATCH is possible */
	static const bool face_match_possible = std::is_same<test_elem_t, trial_elem_t>::value;

	/**
	* \brief return begin iterator of the singular quadrature
	* \return begin iterator of the singular quadrature
	*/
	trial_quadrature_t const &get_trial_quadrature(unsigned idx) const
	{
		return m_face_trial_quadratures[idx];
	}

	/**
	* \brief constructor allocating memory for the quadratures
	*/
	singular_accelerator(void)
	{
		if (face_match_possible)
		{
			for (unsigned idx = 0; idx < test_nset_t::num_nodes; ++idx)
				m_face_trial_quadratures[idx] += trial_blind_t::on_face(
					singular_quadrature_order, test_nset_t::corner_at(idx));
		}
	}

protected:
	/** \brief pointer to the trial singular quadrature */
	trial_quadrature_t m_face_trial_quadratures[test_nset_t::num_nodes];
};


template <
	class TestField, 
	class TrialField, 
	class Singularity, 
	class = typename get_formalism<TestField, TrialField>::type
>
struct double_integral_free_dimensions;

template <class TestField, class TrialField, class Singularity>
struct double_integral_free_dimensions<TestField, TrialField, Singularity, formalism::general>
	: std::integral_constant<int,
	TestField::elem_t::domain_t::dimension + TrialField::elem_t::domain_t::dimension - Singularity::value
> {};


template <class TestField, class TrialField, class Singularity>
struct double_integral_free_dimensions<TestField, TrialField, Singularity, formalism::collocational>
	: std::integral_constant<unsigned,
	TrialField::elem_t::domain_t::dimension
> {};


/** \brief select a singular accelerator for a kernel and test and trial fields
 * \tparam Kernel the kernel type
 * \tparam TestField the test field type
 * \tparam TrialField the trial field tpye
 */
template <class Kernel, class TestField, class TrialField, class Singularity, class = void>
struct select_singular_accelerator
{
	typedef invalid_singular_accelerator type;
};

/** \brief select a singular accelerator for an integrable kernel over a test and trial field
 * \tparam Kernel the kernel type
 * \tparam TestField the test field type
 * \tparam TrialField the trial field tpye
 * \todo Check the enable if relation
 */
template <class Kernel, class TestField, class TrialField, class Singularity>
struct select_singular_accelerator <Kernel, TestField, TrialField, Singularity, typename std::enable_if<
	minimal_reference_dimension<
		typename singular_kernel_traits<Kernel>::singularity_type_t
	>::value <= double_integral_free_dimensions<TestField, TrialField, Singularity>::value
>::type>
{
	typedef singular_accelerator<Kernel, TestField, TrialField> type;
};

}

#endif // SINGULAR_ACCELERATOR_HPP_INCLUDED

