#ifndef PRODUCT_HPP
#define PRODUCT_HPP

template <class A, class B>
struct product_type;

template <> struct product_type<double, double> { typedef double type; };

#endif
