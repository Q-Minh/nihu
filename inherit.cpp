#include <type_traits>
#include <iostream>

class C1
{
public:
	C1(int c1) : m_c1(c1) {}
	int get_c1(void) const { return m_c1; }

private:
	int m_c1;
};

class C2
{
public:
	C2(bool c2) : m_c2(c2) {}
	bool get_c2(void) const { return m_c2; }

private:
	bool m_c2;
};

class C3
{
public:
	C3(double c3) : m_c3(c3) {}
	double get_c3(void) const { return m_c3; }

private:
	double m_c3;
};


/** \brief collect a set of classes by multiple variadic inheritance */
template <class...Args>
class collect : public Args...
{
public:
	/** \brief constructor from individual arguments */
	collect(Args const &...args) : Args(args)...
	{
	}
};


/** \brief merge two collections by preserving unique member types only */
template <class C1, class C2>
struct merge;


/** \brief merge a plain class into a collection if it is already contained */
template <class C1, class D2, bool enable = std::is_base_of<C1, D2>::value>
struct merge_impl
{
	typedef D2 ret_type;

	static ret_type eval(C1 const &c1, D2 const &d2)
	{
		return d2;
	}
};


/** \brief merge a plain class into a collection if it is not contained */
template <class C1, class...Args2>
struct merge_impl<C1, collect<Args2...>, false>
{
	typedef collect<Args2...> d2_type;
	typedef collect<C1, Args2...> ret_type;

	static ret_type eval(C1 const &c1, d2_type const &d2)
	{
		return ret_type(c1, static_cast<Args2 const &>(d2)...);
	}
};


template <class C1, class...Args2>
struct merge<C1, collect<Args2...> > : merge_impl<C1, collect<Args2...> > {};


template <class A1, class...Args2>
struct merge<collect<A1>, collect<Args2...> >
{
	typedef collect<A1> d1_type;
	typedef collect<Args2...> d2_type;
	typedef typename merge<A1, d2_type>::ret_type ret_type;

	static ret_type eval(d1_type const &d1, d2_type const &d2)
	{
		return merge<A1, d2_type>::eval(d1, d2);
	}
};



template <class A1, class...Args1, class...Args2>
struct merge<collect<A1, Args1...>, collect<Args2...> >
{
	typedef collect<A1, Args1...> d1_type;
	typedef collect<Args2...> d2_type;

	typedef typename merge<collect<Args1...>, d2_type>::ret_type inter_type;
	typedef typename merge<A1, inter_type>::ret_type ret_type;

	static ret_type eval(d1_type const &d1, d2_type const &d2)
	{
		return merge<A1, inter_type>::eval(d1,
			merge<collect<Args1...>, d2_type >::eval(d1, d2)
		);
	}
};



int main(void)
{
	std::cout << std::boolalpha;

	C1 c1(-23);
	C2 c2(true);
	C3 c3(2.0);

	collect<C1> d1(c1);
	collect<C2> d2(c2);
	collect<C3> d3(c3);

	collect<C1, C2> d12(c1, c2);
	collect<C2, C3> d23(c2, c3);

	auto d123 = merge<collect<C1>, collect<C2, C3> >::eval(d1, d23);
	auto d123b = merge<collect<C1, C2>, collect<C3> >::eval(d12, d3);
	auto d1223 = merge<collect<C1, C2>, collect<C2, C3> >::eval(d12, d23);

//	auto d1223 = merge(d12, d23);

/*
	collect<C2> d2(c2);
	collect<C1, C2> d12(c1, c2);
	collect<C2, C3> d23(c2, c3);

	auto q = merge(d12, d23);

*/
	std::cout << d1223.get_c1() << d1223.get_c2() << d1223.get_c3();

	return 0;
}

