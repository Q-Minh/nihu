#ifndef HG_HPP_INCLUDED
#define HG_HPP_INCLUDED

#include <iostream>

template <class lhs, class rhs, class op>
struct bin_traits;

struct plus;

template <class T>
struct bin_traits<T, T, plus>
{
	typedef T type;
};

struct plus
{
	template <class lhs, class rhs>
	typename bin_traits<lhs, rhs, plus>::type operator () (lhs const &l, rhs const &r) const
	{
		return l + r;
	}
};

struct index;

template <class Lhs, class Rhs>
struct bin_traits<Lhs, Rhs, index>
{
	typedef typename Lhs::Scalar type;
};

struct index
{
	template <class lhs, class rhs>
	typename bin_traits<lhs, rhs, index>::type operator () (lhs const &l, rhs const &r) const
	{
		return l[r];
	}
};

template <class Derived>
struct hg_traits;

template <class lhs, class rhs, class op>
class hg_bin_expression;

template <class Derived>
class hg_base
{
/*
private:
	hg_base<Derived> &operator=(hg_base<Derived> const &other);
	*/
public:
	typedef typename hg_traits<Derived>::h_t h_t;
	typedef typename hg_traits<Derived>::g_t g_t;

	Derived const &derived(void) const
	{
		return *(static_cast<Derived const *>(this));
	}

	Derived &derived(void)
	{
		return *(static_cast<Derived *>(this));
	}

	template <class otherDerived>
	Derived &operator=(hg_base<otherDerived> const &other)
	{
		derived().get_h() = other.derived().get_h();
		derived().get_g() = other.derived().get_g();
		return derived();
	}

	template <class otherDerived>
	Derived &operator+=(hg_base<otherDerived> const &other)
	{
		derived().get_h() += other.derived().get_h();
		derived().get_g() += other.derived().get_g();
		return derived();
	}

	template <class otherDerived>
	hg_bin_expression<Derived, otherDerived, plus>
		operator +(hg_base<otherDerived> const &other)
	{
		return hg_bin_expression<Derived, otherDerived, plus>(derived(), other.derived(), plus());
	}

	template <class Integer>
	hg_bin_expression<Derived, Integer, index>
		operator [](Integer idx)
	{
		return hg_bin_expression<Derived, Integer, index>(derived(), idx, index());
	}

	std::ostream &print(std::ostream &os) const
	{
		return os << derived().get_h() << '\t' << derived().get_g();
	}
};

template <class Derived>
std::ostream &operator<<(std::ostream &os, hg_base<Derived> const &hg)
{
	return hg.print(os);
}


template <class HType, class GType>
class hg;

template <class HType, class GType>
struct hg_traits<hg<HType, GType> >
{
	typedef HType h_t;
	typedef GType g_t;
};

template <class HType, class GType>
class hg : public hg_base<hg<HType, GType> >
{
public:
	typedef hg_base<hg<HType, GType> > base;
	typedef hg hg_t;
	typedef typename base::h_t h_t;
	typedef typename base::g_t g_t;

	hg(h_t const &h = h_t(), g_t const &g = g_t()) : h(h), g(g)
	{
	}

	// default operator= is explicitly redirected to base's template equivalent
	template <class otherDerived>
	hg &operator=(hg_base<otherDerived> const &other)
	{
		return base::operator=<otherDerived>(other);
	}

	h_t const &get_h(void) const
	{
		return h;
	}

	g_t const &get_g(void) const
	{
		return g;
	}

	h_t &get_h(void)
	{
		return h;
	}

	g_t &get_g(void)
	{
		return g;
	}

protected:
	h_t h;
	g_t g;
};


template <class Lhs, class Rhs, class Op>
struct hg_traits<hg_bin_expression<Lhs, Rhs, Op> >
{
	typedef typename bin_traits<
		typename hg_traits<Lhs>::h_t,
		typename hg_traits<Rhs>::h_t,
		Op
	>::type h_t;
	typedef typename bin_traits<
		typename hg_traits<Lhs>::g_t,
		typename hg_traits<Rhs>::g_t,
		Op
	>::type g_t;
};

template <class Lhs, class Rhs, class Op>
class hg_bin_expression : public hg_base<hg_bin_expression<Lhs, Rhs, Op> >
{
public:
	typedef hg_base<hg_bin_expression<Lhs, Rhs, Op> > base;
	typedef typename base::h_t h_t;
	typedef typename base::g_t g_t;

	hg_bin_expression(Lhs const &lhs, Rhs const &rhs, Op const &op) : lhs(lhs), rhs(rhs), op(op)
	{
	}

	h_t get_h(void) const
	{
		return op(lhs.get_h(), rhs.get_h());
	}

	g_t get_g(void) const
	{
		return op(lhs.get_g(), rhs.get_g());
	}

protected:
	Lhs const &lhs;
	Rhs const &rhs;
	Op op;
};

#endif
