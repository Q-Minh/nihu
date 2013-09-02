/**
 * \file formalism.hpp
 * \brief return weak form formalism from the test and trial field types
 */
#ifndef FORMALISM_HPP
#define FORMALISM_HPP

#include "field.hpp"
#include <type_traits>

/** \brief definition of different weak form formalisms */
namespace formalism
{
	/** \brief general case when the test field is not Dirac */
	struct general {};
	/** \brief collocational case when the test field is Dirac */
	struct collocational {};
}

/** \brief return formalism from Test and Trial field types
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField, class = void>
struct get_formalism;

/** \brief specialiastion of ::get_formalism for the collocational case */
template <class TestField, class TrialField>
struct get_formalism<TestField, TrialField,	typename std::enable_if<
		field_traits<TestField>::is_dirac && !field_traits<TrialField>::is_dirac
	>::type
>
{
	typedef formalism::collocational type;
};


/** \brief specialiastion of ::get_formalism for the general case */
template <class TestField, class TrialField>
struct get_formalism<TestField, TrialField,
	typename std::enable_if<
		!field_traits<TestField>::is_dirac && !field_traits<TrialField>::is_dirac
	>::type
>
{
	typedef formalism::general type;
};


#endif
