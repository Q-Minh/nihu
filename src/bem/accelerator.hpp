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
template <class ElemType>
class accelerator_by
{
public:
	typedef ElemType elem_t;	/**< \brief template argument as nested type */

	typedef typename elem_t::domain_t domain_t;
	typedef gauss_quadrature<domain_t> quadrature_t;
	typedef std::vector<quadrature_t> quadrature_pool_t;

	accelerator_by()
	{
		// initialise quadrature pool
		for (unsigned i = 0; i < 10; ++i)
			m_quadrature_pool.push_back(quadrature_t(i));
	}

	quadrature_pool_t const &get_quadrature_pool(void) const
	{
		return m_quadrature_pool;
	}

protected:
	quadrature_pool_t m_quadrature_pool;
};


template <class ElemTypeVector>
class accelerator
{
public:
	typedef ElemTypeVector elem_type_vector_t;

	template <class elem_t>
	typename accelerator_by<elem_t>::quadrature_pool_t const &
		get_quadrature_pool(void) const
	{
		return container.accelerator_by<elem_t>::get_quadrature_pool();
	}

private:
	template <class T>
	struct acc { typedef accelerator_by<T> type; };

	typename tmp::inherit<
		typename tmp::transform<
		elem_type_vector_t,
		tmp::inserter<tmp::vector<>, tmp::push_back<tmp::_1,tmp::_2> >,
		acc<tmp::_1>
		>::type
	>::type	container;
};

#endif

