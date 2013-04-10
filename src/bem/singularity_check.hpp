#ifndef SINGULARITY_CHECK_HPP_INCLUDED
#define SINGULARITY_CHECK_HPP_INCLUDED

#include <type_traits>
#include "field.hpp"

enum singularity_type {
	REGULAR,
	FACE_MATCH,
	EDGE_MATCH,
	CORNER_MATCH
};

template <class Kernel, class TestField, class TrialField>
class singularity_check
{
public:
	static singularity_type eval(TestField const &test_field, TrialField const &trial_field)
	{
		return REGULAR;
	}
};


template <class Kernel, class Elem, class TrialFieldOption, class TrialDiracOption>
class singularity_check<Kernel, field<Elem, constant_field, dirac_field>, field<Elem, TrialFieldOption, TrialDiracOption> >
{
public:
	typedef field<Elem, constant_field, dirac_field> test_field_t;
	typedef field<Elem, TrialFieldOption, TrialDiracOption>  trial_field_t;

	static bool eval(test_field_t const &test_field, trial_field_t const &trial_field)
	{
		if (test_field.get_elem().get_id() == trial_field.get_elem().get_id())
			return FACE_MATCH;
		return REGULAR;
	}
};


template <class Kernel, class Elem, class TrialField>
class singularity_check<Kernel, field<Elem, constant_field, function_field>, TrialField>
{
public:
	typedef field<Elem, constant_field, function_field> test_field_t;
	typedef TrialField  trial_field_t;
	typedef typename test_field_t::elem_t test_elem_t;
	typedef typename trial_field_t::elem_t trial_elem_t;

	static bool eval(test_field_t const &test_field, trial_field_t const &trial_field)
	{
		// check face match for same element types
		if (std::is_same<test_elem_t, trial_elem_t>::value)
			if (test_field.get_elem().get_id() == trial_field.get_elem().get_id())
				return FACE_MATCH;
		// check edge match

		// check corner match

		return REGULAR;
	}
};





#endif // SINGULARITY_CHECK_HPP_INCLUDED
