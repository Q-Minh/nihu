#ifndef INTEGER_HPP
#define INTEGER_HPP

template <int N>
struct int_
{
	typedef int_ type;

	static int const value = N;
	static int const next = N+1;
	static int const prev = N-1;
};

/**
 * metafunction returning next element
 */
template <class Pos>
struct next;

template <int N>
struct next<int_<N> > : int_<int_<N>::next> {};

/**
 * metafunction returning previous element
 */
template <class Pos>
struct prev;

template <int N>
struct prev<int_<N> > : int_<int_<N>::prev> {};

template <class Pos1, class Pos2>
struct plus;

template <int N, int M>
struct plus<int_<N>, int_<M> > : int_<N+M> {};

template <class Pos1, class Pos2>
struct minus;

template <int N, int M>
struct minus<int_<N>, int_<M> > : int_<N-M> {};

#endif
