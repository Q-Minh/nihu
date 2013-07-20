#ifndef ENABLE_IF_HPP_INCLUDED
#define ENABLE_IF_HPP_INCLUDED

template <class R, bool cond>
struct enable_if { };

template <class R>
struct enable_if<R, true>
{
	typedef R type;
};

#endif // ENABLE_IF_HPP_INCLUDED

