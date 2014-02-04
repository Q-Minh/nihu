#ifndef STORE_OR_ON_THE_FLY_HPP_INCLUDED
#define STORE_OR_ON_THE_FLY_HPP_INCLUDED

#include "core/global_definitions.hpp"

template <bool OnTheFly, class Func, class Ret, class ...Args>
class store_or_on_the_fly;

template <class Func, class Ret, class ...Args>
class store_or_on_the_fly<true, Func, Ret, Args...>
{
public:
	typedef Ret return_type;

	static return_type eval(Args const &...args)
	{
		return Func::eval(args...);
	}
};


template <class Func, class Ret, class...Args>
class store_or_on_the_fly<false, Func, Ret, Args...>
{
public:
	typedef Ret const &return_type;

	static CONSTEXPR return_type eval(Args...)
	{
		return m_stored;
	}

private:
	static const Ret m_stored;
};

template <class Func, class Ret, class ...Args>
const Ret store_or_on_the_fly<false, Func, Ret, Args...>::m_stored = Func::eval(Args()...);

#endif // STORE_OR_ON_THE_FLY_HPP_INCLUDED
