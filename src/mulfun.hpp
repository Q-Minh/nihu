template <int N>
struct int_
{
	typedef int_ type;

	static int const value = N;
	static int const next = N+1;
	static int const prev = N-1;
};

template <class Pos>
struct next;

template <int N>
struct next<int_<N>> : int_<N::next> {};

template <class Pos>
struct prev;

template <int N>
struct prev<int_<N>> : int_<N::prev> {};

template <class Pos1, class Pos2>
struct plus;

template <int N, int M>
struct plus<int_<N>, int_<M>> : int_<N+M> {};

template <class Pos1, class Pos2>
struct minus;

template <int N, int M>
struct minus<int_<N>, int_<M>> : int_<N-M> {};





struct empty{};

struct arg4_tag{};

template <class A0 = empty, class A1 = empty, class A2 = empty, class A3 = empty>
struct arg_seq4
{
	typedef arg_seq4 type; 	/* self-returning */

	typedef arg_seq4_tag tag;	/* tagged */

	typedef A0 arg0;
	typedef A1 arg1;
	typedef A2 arg2;
	typedef A3 arg3;
};


template <class Arg, int Pos>
struct arg_seq4_at;

template <class A0, class A1, class A2, class A3>
struct arg_seq4_at<arg_seq4<A0, A1, A2, A3>, 0> { typedef A0 type; };

template <class A0, class A1, class A2, class A3>
struct arg_seq4_at<arg_seq4<A0, A1, A2, A3>, 1> { typedef A1 type; };

template <class A0, class A1, class A2, class A3>
struct arg_seq4_at<arg_seq4<A0, A1, A2, A3>, 2> { typedef A2 type; };

template <class A0, class A1, class A2, class A3>
struct arg_seq4_at<arg_seq4<A0, A1, A2, A3>, 3> { typedef A3 type; };


template <class tag>
struct at_impl;

template <>
struct at_impl<arg_seq4_tag>
{
	template <class Seq, class Pos>
	struct apply : arg_seq4_at<Seq, Pos::value> {};
};

template <class Seq, class Pos>
struct at : at_impl<typename Seq::tag>::template apply<Seq, Pos> {};


template <class Seq, class Pos>
struct arg_seq4_iterator;

template <class Seq, class Pos>
struct next<iterator<Seq, Pos> > : iterator<Seq, typename next<Pos>::type> {};

template <class Seq, class Pos>
struct prev<iterator<Seq, Pos> > : iterator<Seq, typename prev<Pos>::type> {};

template <class Iter>
struct deref;

template <class Seq, class Pos>
struct deref<Iter<Seq, Pos> > : at<Seq, Pos> {};


