#ifndef QUADRATURE_POOL_HPP_INCLUDED
#define QUADRATURE_POOL_HPP_INCLUDED

#include "kernel.hpp"
#include "field_type_accelerator.hpp"
#include "singular_accelerator.hpp"

template <class Formalism, class Kernel, class TestField, class TrialField>
struct accel_store
{
	typedef singular_accelerator<Formalism, Kernel, TestField, TrialField> singular_accelerator_t;
	static singular_accelerator_t m_singular_accelerator;
};

template<class Formalism, class Kernel, class Test, class Trial>
typename accel_store<Formalism, Kernel, Test, Trial>::singular_accelerator_t
	accel_store<Formalism, Kernel, Test, Trial>::m_singular_accelerator;



template <class Field, class Family>
struct regular_pool_store
{
	typedef field_type_accelerator_pool<Field, Family> pool_t;
	static pool_t m_regular_pool;
};

template <class Field, class Family>
typename regular_pool_store<Field, Family>::pool_t
	regular_pool_store<Field, Family>::m_regular_pool;


#endif // QUADRATURE_POOL_HPP_INCLUDED

