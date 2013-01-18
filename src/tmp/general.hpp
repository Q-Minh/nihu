#ifndef GENERAL_HPP_INCLUDED
#define GENERAL_HPP_INCLUDED

#include "bool.hpp"

template <class A, class B>
struct is_same : false_ {};

template <class A>
struct is_same<A, A> : true_ {};

#endif
