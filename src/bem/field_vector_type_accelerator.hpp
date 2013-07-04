#ifndef FIELD_VECTOR_TYPE_ACCELERATOR_HPP
#define FIELD_VECTOR_TYPE_ACCELERATOR_HPP

#include "field_type_accelerator.hpp"
#include "../tmp/algorithm.hpp"
#include "../tmp/vector.hpp"

template <class FieldVector, class QuadratureFamily>
class field_vector_type_accelerator
{
	typedef FieldVector field_type_vector_t;
	typedef QuadratureFamily quadrature_family_t;

	template <class field_type>
	struct acceleratorize
	{
		typedef field_type_accelerator_pool<field_type, quadrature_family_t> type;
	};

	typedef typename tmp::inherit<
		typename tmp::transform<
		field_type_vector_t,
		tmp::inserter<tmp::vector<>, tmp::push_back<tmp::_1,tmp::_2> >,
		acceleratorize<tmp::_1>
		>::type
	>::type accelerator_container_t;

	template <class field_t>
	typename field_type_accelerator<field_t, quadrature_family_t>::const_iterator
		cbegin(unsigned order) const
	{
		typedef typename acceleratorize<field_t>::type acc_t;
		return m_accelerator_container::acc_t[order]->cbegin();
	}

	template <class field_t>
	typename field_type_accelerator<field_t, quadrature_family_t>::const_iterator
		cend(unsigned order) const
	{
		typedef typename acceleratorize<field_t>::type acc_t;
		return m_accelerator_container::acc_t[order]->cend();
	}

private:
	accelerator_container_t m_accelerator_container;
};

template <class TestFieldVector, class TrialFieldVector, class QuadratureFamily>
class weighted_residual_accelerator
{
public:
	typedef TestFieldVector test_field_vector_t;
	typedef TrialFieldVector trial_field_vector_t;
	typedef QuadratureFamily quadrature_family_t;

	typedef field_vector_type_accelerator<test_field_vector_t, quadrature_family_t> test_accelerator_t;
	typedef field_vector_type_accelerator<trial_field_vector_t, quadrature_family_t> trial_accelerator_t;

private:
	test_accelerator_t m_test_accelerator;
	trial_accelerator_t m_trial_accelerator;
};

#endif // FIELD_VECTOR_TYPE_ACCELERATOR_HPP
