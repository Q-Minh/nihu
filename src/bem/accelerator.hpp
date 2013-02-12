#ifndef ACCELERATOR_HPP_INCLUDED
#define ACCELERATOR_HPP_INCLUDED

#include <vector>

#include "elem_accelerator.hpp"
#include "nset_accelerator.hpp"

#include "../tmp/vector.hpp"
#include "../tmp/control.hpp"

/**
 * \brief an accelerator class storing quadratures and shape sets for an element type and quadrature dimensions
 * \tparam ElemType the accelerated element type
 * \tparam FieldOption field generation option
 * \tparam KernelInput the kernel input class to be stored
 */
template <class ElemType, class FieldOption, class KernelInput>
class accelerator_by
{
public:
	typedef ElemType elem_t;	/**< \brief template argument as nested type */
	typedef FieldOption field_option_t;	/**< \brief template argument as nested type */
	typedef KernelInput kernel_input_t;	/**< \brief template argument as nested type */

	typedef typename elem_t::domain_t domain_t;
	typedef gauss_quadrature<domain_t> quadrature_t;
	typedef typename field<elem_t, field_option_t>::nset_t nset_t;

	typedef std::vector<quadrature_t> quadrature_pool_t;
	typedef nset_pool<nset_t> nset_pool_t;
	typedef elem_pool<elem_t, kernel_input_t> elem_pool_t;

	typedef typename std::vector<elem_pool_t>::const_iterator elem_pool_iterator_t;

	accelerator_by()
	{
		// initialise quadrature pool
		for (unsigned i = 0; i < 10; ++i)
			m_quadrature_pool.push_back(quadrature_t(i));
		// initialise nset-pool
		m_nset_pool = nset_pool_t(m_quadrature_pool.begin(), m_quadrature_pool.end());
	}

	void add_elem(elem_t const &elem)
	{
		elem_pools.push_back(elem_pool_t(elem, m_quadrature_pool.begin(), m_quadrature_pool.end()));
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
		return m_nset_pool;
	}

	quadrature_pool_t const &get_quadrature_pool(void) const
	{
		return m_quadrature_pool;
	}

protected:
	quadrature_pool_t m_quadrature_pool;
	nset_pool_t m_nset_pool;
	std::vector<elem_pool_t> elem_pools;
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
	typename accelerator_by<elem_t, field_option_t, kernel_input_t>::elem_pool_iterator_t
		elem_begin(void) const
	{
		return container.accelerator_by<elem_t, field_option_t, kernel_input_t>::elem_pool_begin();
	}

	template <class elem_t>
	typename accelerator_by<elem_t, field_option_t, kernel_input_t>::elem_pool_iterator_t
		elem_end(void) const
	{
		return container.accelerator_by<elem_t, field_option_t, kernel_input_t>::elem_pool_end();
	}

	template <class elem_t>
	typename accelerator_by<elem_t, field_option_t, kernel_input_t>::nset_pool_t const &
		get_nset_pool(void) const
	{
		return container.accelerator_by<elem_t, field_option_t, kernel_input_t>::get_npool();
	}

	template <class elem_t>
	typename accelerator_by<elem_t, field_option_t, kernel_input_t>::quadrature_pool_t const &
		get_quadrature_pool(void) const
	{
		return container.accelerator_by<elem_t, field_option_t, kernel_input_t>::get_quadrature_pool();
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

