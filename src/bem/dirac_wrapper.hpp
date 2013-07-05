/**
* \file dirac_wrapper.hpp
* \ingroup intop
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief declaration of class ::dirac_wrapper that transforms a ::function_space into a Dirac space
*/
#ifndef DIRAC_WRAPPER_HPP_INCLUDED
#define DIRAC_WRAPPER_HPP_INCLUDED

#include "function_space.hpp"

/**
 *\brief class wrapping a function_space into a Dirac function space
 * \tparam FuncSpace the function space type to wrap
 */
template <class FuncSpace>
class dirac_wrapper
{
public:
	/** \brief template argument as nested type */
	typedef FuncSpace function_space_t;
	
	/**
	 * \brief constructor from a function space reference
	 * \param [in] func_space constant reference to the function space to wrap
	 */
	dirac_wrapper(FuncSpace const &func_space) :
		m_func_space(func_space)
	{
	}
	
	/**
	 * \brief return stored function space reference
	 * \return reference to the wrapped function space
	 */
	function_space_t const &get_function_space(void) const
	{
		return m_func_space;
	}
	
private:
	/** \brief the stored function space reference */
	function_space_t const &m_func_space;
};

/**
 * \brief factory function returning a wrapper fuction space
 * \tparam FuncSpace the type of the function space to wrap
 * \param [in] func_space the function space to wrap
 * \return function space transformed into a Dirac space
 */
template <class FuncSpace>
dirac_wrapper<FuncSpace>
dirac(FuncSpace const &func_space)
{
	return dirac_wrapper<FuncSpace>(func_space);
}

#endif

