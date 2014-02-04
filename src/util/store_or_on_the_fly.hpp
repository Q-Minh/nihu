#ifndef STORE_OR_ON_THE_FLY_HPP_INCLUDED
#define STORE_OR_ON_THE_FLY_HPP_INCLUDED

#include "core/global_definitions.hpp"
#include <type_traits>

template <bool OnTheFly, class Func, class ...Args>
class store_or_on_the_fly;

template <class Func, class ...Args>
class store_or_on_the_fly<true, Func, Args...>
{
public:
	typedef decltype(Func::eval(Args()...)) functor_ret_type;
	typedef functor_ret_type return_type;

	static return_type eval(Args const &...args)
	{
		return Func::eval(args...);
	}
};


template <class Func, class...Args>
class store_or_on_the_fly<false, Func, Args...>
{
public:
	typedef decltype(Func::eval(Args()...)) functor_ret_type;
	typedef functor_ret_type const &return_type;

	static CONSTEXPR return_type eval(Args...)
	{
		return m_stored;
	}

private:
	static const functor_ret_type m_stored;
};

template <class Func, class ...Args>
const typename store_or_on_the_fly<false, Func, Args...>::functor_ret_type
store_or_on_the_fly<false, Func, Args...>::m_stored = Func::eval(Args()...);

#endif // STORE_OR_ON_THE_FLY_HPP_INCLUDED
