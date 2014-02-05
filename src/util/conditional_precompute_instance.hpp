#ifndef CONDITIONAL_PRECOMPUTE_INSTANCE_HPP_INCLUDED
#define CONDITIONAL_PRECOMPUTE_INSTANCE_HPP_INCLUDED

#include "core/global_definitions.hpp"
#include <type_traits>

/** \brief Conditionally precompute and store objects
 * \tparam OnThFly true if the object does not need to be precomputed
 * \tparam Func the class that computes object instances
 * \tparam Arguments passed to the Func class
 */
template <bool OnTheFly, class Func, class ...Args>
class conditional_precompute_instance
{
public:
	/** \brief the return type of the Func function class */
	typedef decltype(Func::eval(Args()...)) functor_ret_type;
	/** \brief the return type of the static function eval */
	typedef functor_ret_type return_type;

	template <class...ConstrArgs>
	conditional_precompute_instance(ConstrArgs const &...)
	{
	}

	/** \brief return object computed by the function class */
	return_type operator()(Args const &...args) const
	{
		return Func::eval(args...);
	}
};


template <class Func, class...Args>
class conditional_precompute_instance<false, Func, Args...>
{
public:
	/** \brief the return type of the Func function class */
	typedef decltype(Func::eval(Args()...)) functor_ret_type;
	/** \brief the return type of the static function eval */
	typedef typename std::add_lvalue_reference<
        typename std::add_const<functor_ret_type>::type
    >::type return_type;

	conditional_precompute_instance(Args const &...args)
        : m_stored(Func::eval(args...))
	{
	}

	/** \brief return object computed by the function class */
	template <class...EvalArgs>
	return_type operator()(EvalArgs const &...) const
	{
		return m_stored;
	}

private:
    const functor_ret_type m_stored;
};


#endif // CONDITIONAL_PRECOMPUTE_INSTANCE_HPP_INCLUDED
