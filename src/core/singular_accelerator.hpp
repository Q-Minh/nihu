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
 * \file singular_accelerator.hpp
 * \ingroup quadrature
 * \brief implementation of class ::singular_accelerator
 */

#ifndef SINGULAR_ACCELERATOR_HPP_INCLUDED
#define SINGULAR_ACCELERATOR_HPP_INCLUDED

#include "kernel.hpp"	// for the galerkin case
#include "singular_galerkin_quadrature.hpp"	// for the galerkin case
#include "blind_transform_selector.hpp"	// for the collocation case
#include "blind_singular_quadrature.hpp"	// for the collocation case
#include "field.hpp"

#include "../util/dual_range.hpp"
#include "formalism.hpp"

#include <stdexcept>

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
template <class Kernel, class TestField, class TrialField, class = typename get_formalism<TestField, TrialField>::type >
class singular_accelerator;


/**
* \brief specialisation of ::singular_accelerator for the general formalism
* \tparam Kernel the kernel that is integrated
* \tparam TestField the test field over which integration is performed
* \tparam TrialField the trial field over which integration is performed
*/
template <class Kernel, class TestField, class TrialField>
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

	/** \brief the L-set of the test domain */
	typedef typename test_field_t::elem_t::lset_t test_lset_t;
	/** \brief the L-set of the trial domain */
	typedef typename trial_field_t::elem_t::lset_t trial_lset_t;

	/** \brief test domain */
	typedef typename test_field_t::elem_t::domain_t test_domain_t;
	/** \brief trial domain */
	typedef typename trial_field_t::elem_t::domain_t trial_domain_t;

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

	/** \brief indicates whether 2d match is possible */
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
		case 2:
			return iterator(
				m_face_test_quadrature.begin(),
				m_face_trial_quadrature.begin());
			break;
		case 1:
		{
			unsigned test_domain_corner =
				test_lset_t::node_to_domain_corner(elem_match.get_overlap().get_ind1());
			unsigned trial_domain_corner =
				trial_lset_t::node_to_domain_corner(elem_match.get_overlap().get_ind2());
			return iterator(
				m_edge_test_quadrature[test_domain_corner].begin(),
				m_edge_trial_quadrature[trial_domain_corner].begin());
			break;
		}
		case 0:
		{
			unsigned test_domain_corner =
				test_lset_t::node_to_domain_corner(elem_match.get_overlap().get_ind1());
			unsigned trial_domain_corner =
				trial_lset_t::node_to_domain_corner(elem_match.get_overlap().get_ind2());
			return iterator(
				m_corner_test_quadrature[test_domain_corner].begin(),
				m_corner_trial_quadrature[trial_domain_corner].begin());
			break;
		}
		case -1:
			throw std::logic_error("Cannot return singular quadrature for regular type");
			break;
		default:
			throw std::invalid_argument("Unknown singularity type");
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
		case 2:
			return iterator(
				m_face_test_quadrature.end(),
				m_face_trial_quadrature.end());
			break;
		case 1:
		{
			unsigned test_domain_corner =
				test_lset_t::node_to_domain_corner(elem_match.get_overlap().get_ind1());
			unsigned trial_domain_corner =
				trial_lset_t::node_to_domain_corner(elem_match.get_overlap().get_ind2());
			return iterator(
				m_edge_test_quadrature[test_domain_corner].end(),
				m_edge_trial_quadrature[trial_domain_corner].end());
			break;
		}
		case 0:
		{
			unsigned test_domain_corner =
				test_lset_t::node_to_domain_corner(elem_match.get_overlap().get_ind1());
			unsigned trial_domain_corner =
				trial_lset_t::node_to_domain_corner(elem_match.get_overlap().get_ind2());
			return iterator(
				m_corner_test_quadrature[test_domain_corner].end(),
				m_corner_trial_quadrature[trial_domain_corner].end());
			break;
		}
		case -1:
			throw std::logic_error("Cannot return singular quadrature for NO_MATCH type");
			break;
		default:
			throw std::invalid_argument("Unknown singularity type");
		}
	}

private:
	void generate_face(std::false_type)
	{
	}

	void generate_face(std::true_type)
	{
		// construct facials
		quad_factory_t::template generate<match::match_2d_type>(
			m_face_test_quadrature,
			m_face_trial_quadrature,
			singular_quadrature_order);
	}

public:
	/**
	* \brief constructor allocating space for the quadratures
	*/
	singular_accelerator(void)
	{
		// generate face quadratures with separate function
		// so that it does not compile of not needed (different domains)
		generate_face(std::integral_constant<bool, face_match_possible>());

		// create temporary quadratures for rotating
		test_quadrature_t test_edge_q, test_corner_q;
		trial_quadrature_t trial_edge_q, trial_corner_q;

		// create test quadrature singular on the first corner
		quad_factory_t::template generate<match::match_0d_type>(
			test_corner_q, trial_corner_q, singular_quadrature_order);
		// create test quadrature singular on the first edge
		quad_factory_t::template generate<match::match_1d_type>(
			test_edge_q, trial_edge_q, singular_quadrature_order);

		// rotate test quadratures
		for	(unsigned i = 0; i < n_test_corners; ++i)
		{
			// fill transform coordinates
			Eigen::Matrix<scalar_t, n_test_corners, test_domain_t::dimension> test_corners;
			for (unsigned k = 0; k < n_test_corners; ++k)
				test_corners.row(k) = test_domain_t::get_corner((i+k) % n_test_corners);

			// rotate
			m_corner_test_quadrature[i] =
				test_corner_q.template transform<isoparam_shape_set<test_domain_t> >(test_corners);
			m_edge_test_quadrature[i] =
				test_edge_q.template transform<isoparam_shape_set<test_domain_t> >(test_corners);
		}
		// rotate trial quads
		for (unsigned i = 0; i < n_trial_corners; ++i)
		{
			// fill transform coordinates
			Eigen::Matrix<scalar_t, n_trial_corners, trial_domain_t::dimension> trial_corners;
			for (unsigned k = 0; k < n_trial_corners; ++k)
				trial_corners.row(k) = trial_domain_t::get_corner((i+k) % n_trial_corners);

			// rotate
			m_corner_trial_quadrature[i] =
				trial_corner_q.template transform<isoparam_shape_set<trial_domain_t> >(trial_corners);

			// when dealing with the EDGE_MATCH case, we need to take into consideration that the test and
			// trial elements contain the singular edge in opposite nodal order. Therefore, if the
			// test element's singular edge is 0-1, then the trial element's singular edge should be 1-0
			// this numbering is implemented below: i+1-k instead of i+k

			// fill transform coordinates
			for (unsigned k = 0; k < n_trial_corners; ++k)
				trial_corners.row(k) =
				trial_domain_t::get_corner((i+1+n_trial_corners-k) % n_trial_corners);

			// rotate but keep weights unadjusted
			m_edge_trial_quadrature[i] = trial_edge_q.template transform<isoparam_shape_set<trial_domain_t> >(trial_corners);
			m_edge_trial_quadrature[i] *= -1.0;
		}
	}

protected:
	/**\brief singular quadrature on the test elem for FACE_MATCH case */
	test_quadrature_t m_face_test_quadrature;
	/**\brief singular quadrature on the trial elem for FACE_MATCH case */
	trial_quadrature_t m_face_trial_quadrature;

	/**\brief singular quadratures on the test elem for CORNER_MATCH case */
	test_quadrature_t m_corner_test_quadrature[test_domain_t::num_corners];
	/**\brief singular quadratures on the trial elem for CORNER_MATCH case */
	trial_quadrature_t m_corner_trial_quadrature[trial_domain_t::num_corners];

	/**\brief singular quadratures on the test elem for EDGE_MATCH case */
	test_quadrature_t m_edge_test_quadrature[test_domain_t::num_corners];
	/**\brief singular quadratures on the trial elem for EDGE_MATCH case */
	trial_quadrature_t m_edge_trial_quadrature[trial_domain_t::num_corners];
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


/** \brief select a singular accelerator for a kernel and test and trial fields
 * \tparam Kernel the kernel type
 * \tparam TestField the test field type
 * \tparam TrialField the trial field tpye
 */
template <class Kernel, class TestField, class TrialField, class = void>
struct select_singular_accelerator
{
	typedef invalid_singular_accelerator type;
};

/** \brief select a singular accelerator for an integrable kernel over a test and trial field
 * \tparam Kernel the kernel type
 * \tparam TestField the test field type
 * \tparam TrialField the trial field tpye
 */
template <class Kernel, class TestField, class TrialField>
struct select_singular_accelerator <Kernel, TestField, TrialField, typename std::enable_if<
	minimal_reference_dimension<
		typename singular_kernel_traits<Kernel>::singularity_type_t
	>::value <= TrialField::elem_t::domain_t::dimension
>::type>
{
	typedef singular_accelerator<Kernel, TestField, TrialField> type;
};

#endif // SINGULAR_ACCELERATOR_HPP_INCLUDED

