/** \file match_types.hpp */

#ifndef MATCH_TYPES_HPP_INCLUDED
#define MATCH_TYPES_HPP_INCLUDED

#include <type_traits>
#include "formalism.hpp"
#include "element_match.hpp"
#include "../tmp/algorithm.hpp"

/** \brief matafunction assigning a match type vector to two fields
 * \details The match type vector is a compile time vector containing
 * possible match dimensions.
 */
template <class TestField, class TrialField, class Enable = void>
struct match_type_vector;

/** \brief specialisation of match_type_vector for the general formalism
 *  \todo incorporate that a d-match is only possible if the element domains are the same
 */
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

public:
	typedef typename std::conditional<
		d == 2,
		tmp::vector<match::match_0d_type, match::match_1d_type, match::match_2d_type>,
		typename std::conditional<
			d == 1,
			tmp::vector<match::match_0d_type, match::match_1d_type>,
			typename std::conditional<
				d == 0,
				tmp::vector<match::match_0d_type>,
				tmp::vector<>
			>::type
		>::type
	>::type type;
};

/** \brief specialisation of match_type_vector for the collocational formalism
 * \todo incorporate that a d-match is only possible if the element domains are the same
 */
template <class TestField, class TrialField>
struct match_type_vector<TestField, TrialField,
	typename std::enable_if<
		std::is_same<
			typename get_formalism<TestField, TrialField>::type,
			formalism::collocational
		>::value
	>::type
>
{
private:
	enum { d = TestField::elem_t::domain_t::dimension };

	// start with empty vector
	typedef tmp::vector<> t;

	// push back 0d if needed
	typedef typename std::conditional<
		tmp::is_member<typename shape_set_traits::position_dof_vector<typename TestField::nset_t>::type, dof0>::value
		&&
		tmp::is_member<typename shape_set_traits::position_dof_vector<typename TrialField::nset_t>::type, dof0>::value,
		typename tmp::push_back<t, match::match_0d_type>::type,
		t
	>::type type0;

	// push back 0d if needed
	typedef typename std::conditional<
		tmp::is_member<typename shape_set_traits::position_dof_vector<typename TestField::nset_t>::type, dof1>::value
		&&
		tmp::is_member<typename shape_set_traits::position_dof_vector<typename TrialField::nset_t>::type, dof1>::value,
		typename tmp::push_back<type0, match::match_1d_type>::type,
		type0
	>::type type1;

	// push back 2d if needed
	typedef typename std::conditional<
		tmp::is_member<typename shape_set_traits::position_dof_vector<typename TestField::nset_t>::type, dof2>::value
		&&
		tmp::is_member<typename shape_set_traits::position_dof_vector<typename TrialField::nset_t>::type, dof2>::value,
		typename tmp::push_back<type1, match::match_2d_type>::type,
		type1
	>::type type2;

public:
	typedef type2 type;
};

#endif // MATCH_TYPES_HPP_INCLUDED

