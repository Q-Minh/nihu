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

/** \file match_types.hpp */

#ifndef MATCH_TYPES_HPP_INCLUDED
#define MATCH_TYPES_HPP_INCLUDED

#include <type_traits>
#include "formalism.hpp"
#include "../tmp/algorithm.hpp"

/** \brief match types */
namespace match
{
	/** \brief no singularity */
	typedef std::integral_constant<int, -1> no_match_type;

	/** \brief two elements are adjacent with a 0d (vertex) match */
	typedef std::integral_constant<int, 0> match_0d_type;
	/** \brief two elements are adjacent with a 1d (line) match */
	typedef std::integral_constant<int, 1> match_1d_type;
	/** \brief two elements are adjacent with a 2d (surface) match */
	typedef std::integral_constant<int, 2> match_2d_type;
}

/** \brief matafunction assigning a match type vector to two fields
 * \details The match type vector is a compile time vector containing
 * possible match dimensions.
 */
template <class TestField, class TrialField, class Enable = void>
struct match_type_vector;

/** \brief specialisation of match_type_vector for the general formalism */
template <class TestField, class TrialField>
struct match_type_vector<TestField, TrialField,
	typename std::enable_if<
		std::is_same<
			typename get_formalism<TestField, TrialField>::type,
			formalism::general
		>::value
	>::type
>
{
private:
	enum { d = TestField::elem_t::domain_t::dimension };
	typedef typename TestField::elem_t test_elem_t;
	typedef typename TrialField::elem_t trial_elem_t;

	template <class vector, unsigned dim>
	struct push_if_needed : std::conditional<
		(d > dim) ||
		(d == dim && std::is_same<test_elem_t, trial_elem_t>::value),
		typename tmp::push_back<vector, std::integral_constant<int, int(dim)> >::type,
		vector
	> {};

public:
	typedef typename push_if_needed<
		typename push_if_needed<
			typename push_if_needed<
				tmp::vector<>, 0>::type, 1
		>::type, 2
	>::type type;
};

/** \brief specialisation of match_type_vector for the collocational formalism */
template <class TestField, class TrialField>
struct match_type_vector<TestField, TrialField, typename std::enable_if<
	std::is_same<
		typename get_formalism<TestField, TrialField>::type,
		formalism::collocational
	>::value
>::type>
{
private:
	enum { d = TestField::elem_t::domain_t::dimension };
	typedef typename TestField::elem_t test_elem_t;
	typedef typename shape_set_traits::position_dof_vector<typename TestField::nset_t>::type test_dof_vector;
	typedef typename TrialField::elem_t trial_elem_t;

	template <class vector, unsigned dim>
	struct push_if_needed : std::conditional<
		(d > dim && tmp::is_member<test_dof_vector, position_dof<dim> >::value) ||
		(d == dim && std::is_same<test_elem_t, trial_elem_t>::value &&
		 tmp::is_member<test_dof_vector, position_dof<dim> >::value),
		typename tmp::push_back<vector, std::integral_constant<int, int(dim)> >::type,
		vector
	> {};

public:
	typedef typename push_if_needed<
		typename push_if_needed<
			typename push_if_needed<
				tmp::vector<>, 0>::type, 1
		>::type, 2
	>::type type;
};

#endif // MATCH_TYPES_HPP_INCLUDED

