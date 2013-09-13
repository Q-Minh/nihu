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
 * \file element_match.hpp
 * \ingroup quadrature
 * \brief determine element singularities
 */

#ifndef ELEMENT_MATCH_HPP_INCLUDED
#define ELEMENT_MATCH_HPP_INCLUDED

#include "element.hpp"

/** \brief singularity types */
enum singularity_type {
	/** \brief no singularity */
	REGULAR,
	/** \brief two elements are identical */
	FACE_MATCH,
	/** \brief two elements share a common edge */
	EDGE_MATCH,
	/** \brief two elements share a common corner */
	CORNER_MATCH
};

/** \brief singularity types */
namespace singularity
{
	/** \brief no singularity */
	typedef std::integral_constant<unsigned, REGULAR> regular_type;
	/** \brief two elements are identical */
	typedef std::integral_constant<unsigned, REGULAR> face_match_type;
	/** \brief two elements share a common edge */
	typedef std::integral_constant<unsigned, REGULAR> edge_match_type;
	/** \brief two elements share a common corner */
	typedef std::integral_constant<unsigned, REGULAR> corner_match_type;
}

/** \brief class describing the adjacency (match) state of two elements */
class element_match
{
public:
	/** \brief constructor
	 * \param [in] sing_type the singularity type
	 * \param [in] overlap the overlapping state
	 */
	element_match(
		singularity_type const &sing_type,
		element_overlapping const &overlap = element_overlapping()) :
		m_sing_type(sing_type), m_overlap(overlap)
	{
	}

	/** \brief return singularity type */
	singularity_type const &get_singularity_type(void) const
	{
		return m_sing_type;
	}

	/** \brief return overlapping state */
	element_overlapping const &get_overlap(void) const
	{
		return m_overlap;
	}

private:
	/** \brief the singularity type */
	singularity_type const m_sing_type;
	/** \brief the overlapping state */
	element_overlapping const m_overlap;
};


/** \brief function to determine the overlapping state of two elements
 * \tparam test_field_t the test field type
 * \tparam trial_field_t the trial field type
 * \param [in] test_field the test field
 * \param [in] trial_field the trial field
 * \return the element match object
 */
template <class test_field_t, class trial_field_t>
element_match element_match_eval(
	field_base<test_field_t> const &test_field,
	field_base<trial_field_t> const &trial_field)
{
	bool const face_match_possible = std::is_same<
		typename test_field_t::elem_t,
		typename trial_field_t::elem_t
	>::value;

	if (face_match_possible) // compile time if
		if (test_field.get_elem().get_id() == trial_field.get_elem().get_id())
			return element_match(FACE_MATCH);

	// compile time if, only for the general approach
	if (std::is_same<typename get_formalism<test_field_t, trial_field_t>::type, formalism::general>::value)
	{
		element_overlapping overlap(test_field.get_elem().get_overlapping(trial_field.get_elem()));

		if (overlap.get_num() > 1)
			return element_match(EDGE_MATCH, overlap);

		if (overlap.get_num() == 1)
			return element_match(CORNER_MATCH, overlap);
	}

	return element_match(REGULAR);
}

#endif // ELEMENT_MATCH_HPP_INCLUDED
