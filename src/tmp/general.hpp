#ifndef GENERAL_HPP_INCLUDED
#define GENERAL_HPP_INCLUDED

#include "bool.hpp"

template <class A, class B>
struct is_same : false_ {};

template <class A>
struct is_same<A, A> : true_ {};


template <class T>
struct is_reference : false_ {};

template <class T>
struct is_reference<T &> : true_ {};

template <class T>
struct add_ref
{
	typedef typename if_<
		typename is_reference<T>::type,
		T,
		T &
	>::type type;
};


#endif
