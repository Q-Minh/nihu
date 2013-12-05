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
#include "field.hpp"
#include "formalism.hpp"

/** \brief element match type */
enum element_match_type {
	/** \brief no singularity */
	NO_MATCH,
	/** \brief two elements are identical */
	FACE_MATCH,
	/** \brief two elements share a common edge */
	EDGE_MATCH,
	/** \brief two elements share a common corner */
	CORNER_MATCH
};

/** \brief match types */
namespace match
{
	/** \brief no singularity */
	typedef std::integral_constant<unsigned, NO_MATCH> no_match_type;
	/** \brief two elements are identical */
	typedef std::integral_constant<unsigned, FACE_MATCH> face_match_type;
	/** \brief two elements share a common edge */
	typedef std::integral_constant<unsigned, EDGE_MATCH> edge_match_type;
	/** \brief two elements share a common corner */
	typedef std::integral_constant<unsigned, CORNER_MATCH> corner_match_type;
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
		element_match_type const &match_type,
		element_overlapping const &overlap = element_overlapping()) :
		m_match_type(match_type), m_overlap(overlap)
	{
	}

	/** \brief return singularity type */
	element_match_type const &get_singularity_type(void) const
	{
		return m_match_type;
	}

	/** \brief return overlapping state */
	element_overlapping const &get_overlap(void) const
	{
		return m_overlap;
	}

private:
	/** \brief the singularity type */
	element_match_type const m_match_type;
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

	return element_match(NO_MATCH);
}

#endif // ELEMENT_MATCH_HPP_INCLUDED
