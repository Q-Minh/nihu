#ifndef ACCELERATOR_HPP_INCLUDED
#define ACCELERATOR_HPP_INCLUDED

#include <vector>

#include "elem_accelerator.hpp"
#include "nset_accelerator.hpp"

#include "../tmp/vector.hpp"
#include "../tmp/control.hpp"

template <class ElemType, class FieldOption, class KernelInput>
class accelerator_by
{
public:
	typedef ElemType elem_t;
	typedef FieldOption field_option_t;
	typedef KernelInput kernel_input_t;

	typedef typename elem_t::domain_t domain_t;
	typedef gauss_quadrature<domain_t> quadrature_t;
	typedef typename field<elem_t, field_option_t>::nset_t nset_t;

	typedef elem_pool<elem_t, kernel_input_t> elem_pool_t;
	typedef nset_pool<nset_t> nset_pool_t;

	typedef typename std::vector<elem_pool_t>::const_iterator elem_pool_iterator_t;

	accelerator_by()
	{
		// initialise quadrature pool
		q.push_back(quadrature_t(1));
		q.push_back(quadrature_t(2));
		q.push_back(quadrature_t(3));
		// initialise nset-pool
		n_pool = nset_pool_t(q.begin(), q.end());
	}

	void add_elem(elem_t const &elem)
	{
		elem_pools.push_back(elem_pool_t(elem, q.begin(), q.end()));
	}

	elem_pool_iterator_t elem_pool_begin(void) const
	{
		return elem_pools.begin();
	}

	elem_pool_iterator_t elem_pool_end(void) const
	{
		return elem_pools.end();
	}

	nset_pool_t const &get_npool(void) const
	{
		return n_pool;
	}

protected:
	std::vector<quadrature_t> q;
	std::vector<elem_pool_t> elem_pools;
	nset_pool_t n_pool;
};


template <class ElemTypeVector, class FieldOption, class KernelInput>
class accelerator
{
public:
	typedef ElemTypeVector elem_type_vector_t;
	typedef FieldOption field_option_t;
	typedef KernelInput kernel_input_t;

	template <class elem_t>
	void add_elem(elem_t const &elem)
	{
		container.accelerator_by<elem_t, field_option_t, kernel_input_t>::add_elem(elem);
	}

	template <class elem_t>
	struct nset_pool_t { typedef typename accelerator_by<elem_t, field_option_t, kernel_input_t>::nset_pool_t type; };

	template <class elem_t>
	struct elem_pool_t { typedef typename accelerator_by<elem_t, field_option_t, kernel_input_t>::elem_pool_t type; };

	template <class elem_t>
	typename accelerator_by<elem_t, field_option_t, kernel_input_t>::elem_pool_iterator_t elem_begin(void) const
	{
		return container.accelerator_by<elem_t, field_option_t, kernel_input_t>::elem_pool_begin();
	}

	template <class elem_t>
	typename accelerator_by<elem_t, field_option_t, kernel_input_t>::elem_pool_iterator_t elem_end(void) const
	{
		return container.accelerator_by<elem_t, field_option_t, kernel_input_t>::elem_pool_end();
	}

	template <class elem_t>
	typename accelerator_by<elem_t, field_option_t, kernel_input_t>::nset_pool_t const &
	get_nset_pool(void) const
	{
		return container.accelerator_by<elem_t, field_option_t, kernel_input_t>::get_npool();
	}

private:
	template <class T>
	struct acc { typedef accelerator_by<T, field_option_t, kernel_input_t> type; };

	typename tmp::inherit<
		typename tmp::transform<
		elem_type_vector_t,
		tmp::inserter<tmp::vector<>, tmp::push_back<tmp::_1,tmp::_2> >,
		acc<tmp::_1>
		>::type
	>::type	container;
};

#endif

