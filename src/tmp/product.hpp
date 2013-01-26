#ifndef PRODUCT_HPP
#define PRODUCT_HPP

template <class A, class B>
struct product_type;

template <class T>
struct product_type<T, T> { typedef T type; };

#endif

