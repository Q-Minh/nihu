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
#include <type_traits>

/** \brief match types */
namespace match
{
	/** \brief no singularity */
	typedef std::integral_constant<int, -1> no_match_type;
	/** \brief two elements share a common corner */
	typedef std::integral_constant<int, 0> match_0d_type;
	/** \brief two elements share a common edge */
	typedef std::integral_constant<int, 1> match_1d_type;
	/** \brief two elements are identical */
	typedef std::integral_constant<int, 2> match_2d_type;
}

/** \brief class describing the adjacency (match) state of two elements */
class element_match
{
public:
	/** \brief constructor
	 * \param [in] match_dimension the match dimension
	 * \param [in] overlap the overlapping state
	 */
	element_match(
		int const match_dimension,
		element_overlapping const &overlap = element_overlapping()) :
		m_match_dimension(match_dimension), m_overlap(overlap)
	{
	}

	/** \brief return singularity type */
	int get_match_dimension(void) const
	{
		return m_match_dimension;
	}

	/** \brief return overlapping state */
	element_overlapping const &get_overlap(void) const
	{
		return m_overlap;
	}

private:
	/** \brief the singularity type */
	int m_match_dimension;
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
			return element_match(2);

	// compile time if, only for the general approach
	if (std::is_same<typename get_formalism<test_field_t, trial_field_t>::type, formalism::general>::value)
	{
		element_overlapping overlap(test_field.get_elem().get_overlapping(trial_field.get_elem()));

		if (overlap.get_num() > 1)
			return element_match(1, overlap);

		if (overlap.get_num() == 1)
			return element_match(0, overlap);
	}

	// no match
	return element_match(-1);
}

#endif // ELEMENT_MATCH_HPP_INCLUDED

