#ifndef FORMALISM_HPP
#define FORMALISM_HPP

#include "field.hpp"
#include <type_traits>

namespace formalism
{
	struct general {};
	struct collocational {};
}

template <class TestField, class TrialField, class = void>
struct get_formalism;

template <class TestField, class TrialField>
struct get_formalism<TestField, TrialField,	typename std::enable_if<
		field_traits<TestField>::is_dirac && !field_traits<TrialField>::is_dirac
	>::type
>
{
	typedef formalism::collocational type;
};


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
