#ifndef BOOL_HPP
#define BOOL_HPP

template <bool x>
struct bool_
{
	static bool const value = x;
	typedef bool_ type;
};

typedef bool_<true> true_;
typedef bool_<false> false_;


template <class C, class T, class F>
struct if_;

template <class T, class F>
struct if_<true_, T, F> { typedef T type; };

template <class T, class F>
struct if_<false_, T, F> { typedef F type; };


template <class A> struct not_ : bool_<!A::value> {};
template <class A, class B, class C = false_> struct or_ : bool_<A::value || B::value || C::value> {};
template <class A, class B, class C = true_> struct and_ : bool_<A::value && B::value && C::value> {};

#endif

