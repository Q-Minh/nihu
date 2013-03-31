/**
 * \file accelerator.hpp
 * \brief Declaration of class accelerator
 * \author Peter Fiala and Peter Rucz, fiala@hit.bme.hu, rucz@hit.bme.hu
 */

#ifndef ACCELERATOR_HPP_INCLUDED
#define ACCELERATOR_HPP_INCLUDED

#include <vector>

#include "../tmp/vector.hpp"
#include "../tmp/control.hpp"

/**
 * \brief an accelerator class storing quadratures for an element type
 * \tparam ElemType the accelerated element type
 */
template <class ElemType>
class accelerator_by
{
public:
	typedef ElemType elem_t;	/**< \brief template argument as nested type */

	static const size_t MAX_ORDER = 10;	/**< \brief maximal quadrature order */

	typedef typename elem_t::domain_t domain_t;	/**< \brief the domain type of the element */
	typedef gauss_quadrature<domain_t> quadrature_t;	/**< \brief the quadrature type of the element */
	typedef std::vector<quadrature_t> quadrature_pool_t;	/**< \brief the quadrature pool type of the element */

	/** \brief constructor initialising the quadratures */
	accelerator_by()
	{
		// initialise quadrature pool
		for (unsigned i = 0; i < MAX_ORDER; ++i)
			m_quadrature_pool.push_back(quadrature_t(i));
	}

	/** \brief return reference to the quadrature pool
	 * \return reference to the quadrature pool
	 */
	quadrature_pool_t const &get_quadrature_pool(void) const
	{
		return m_quadrature_pool;
	}

protected:
	quadrature_pool_t m_quadrature_pool;	/**< \brief all quadratures of the element type up to order MAX_ORDER */
};


/**
 * \brief an accelerator class storing quadratures for a vector of element types
 * \tparam ElemTypeVector vector of element types to accelerate
 */
template <class ElemTypeVector>
class accelerator
{
public:
	typedef ElemTypeVector elem_type_vector_t; /**< \brief template argument as nested type */

	/** \brief return quadrature pool of a specific element type
	 * \tparam elem_t the element type
	 * \return quadrature pool of a specific element type
	 */
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

#endif // ifndef ACCELERATOR_HPP_INCLUDED

