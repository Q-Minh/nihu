#ifndef ITERATOR_HPP
#define ITERATOR_HPP

template <class Iter>
struct iterator_value_type
{
	typedef typename Iter::value_type type;
};

template <class T>
struct iterator_value_type<T*>
{
	typedef T type;
};

#endif

