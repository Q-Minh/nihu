#ifndef CONDITIONAL_PRECOMPUTE_HPP_INCLUDED
#define CONDITIONAL_PRECOMPUTE_HPP_INCLUDED

#include "core/global_definitions.hpp"
#include <type_traits>

/** \brief Conditionally precompute and store objects
 * \tparam OnThFly true if the object does not need to be precomputed
 * \tparam Func the class that computes object instances
 * \tparam Arguments passed to the Func class
 */
template <bool OnTheFly, class Func, class ...Args>
class conditional_precompute
{
public:
	/** \brief the return type of the Func function class */
	typedef decltype(Func::eval(Args()...)) functor_ret_type;
	/** \brief the return type of the static function eval */
	typedef functor_ret_type return_type;

	/** \brief return object computed by the function class */
	static return_type eval(Args const &...args)
	{
		return Func::eval(args...);
	}
};


template <class Func, class...Args>
class conditional_precompute<false, Func, Args...>
{
public:
	/** \brief the return type of the Func function class */
	typedef decltype(Func::eval(Args()...)) functor_ret_type;
	/** \brief the return type of the static function eval */
	typedef functor_ret_type const &return_type;

	/** \brief return object computed by the function class */
	static return_type eval(Args const &...)
	{
		static const functor_ret_type m_stored = Func::eval(Args()...);
		return m_stored;
	}
};


#endif // CONDITIONAL_PRECOMPUTE_HPP_INCLUDED
