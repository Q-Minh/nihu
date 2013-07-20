#ifndef NEWCOUPLE_HPP_INCLUDED
#define NEWCOUPLE_HPP_INCLUDED

#include "enable_if.hpp"
#include "type_info.hpp"
#include <type_traits>
#include <utility> // pair
#include <iostream>

template <class T, class S = T>
class couple :
	private std::pair<T, S>
{
public:
	typedef std::true_type is_couple;

	typedef std::pair<T, S> base_t;

	couple(T &&left, S &&right) :
		base_t(std::forward<T>(left), std::forward<S>(right))
	{
		std::cout << "constructing couple< ";
		print_type_info<T>();
		std::cout << " , ";
		print_type_info<S>();
		std::cout << ">\n";
	}

	T &get_first()
	{
		return base_t::first;
	}

	T const &get_first() const
	{
		return base_t::first;
	}

	S &get_second()
	{
		return base_t::second;
	}

	S const &get_second() const
	{
		return base_t::second;
	}
};


template <class LC, class RC>
class couple_product
{
private:
	LC m_lhs;
	RC m_rhs;

public:
	couple_product(LC &&lhs, RC &&rhs) :
		m_lhs(std::forward<LC>(lhs)),
		m_rhs(std::forward<RC>(rhs))
	{
	}

	auto get_first(void) const -> decltype(m_lhs.get_first() * m_rhs.get_first())
	{
		return m_lhs.get_first() * m_rhs.get_first();
	}

	auto get_second(void) const -> decltype(m_lhs.get_second() * m_rhs.get_second())
	{
		return m_lhs.get_second() * m_rhs.get_second();
	}
};


template <class C, class Rhs>
class couple_product_right
{
private:
	C m_couple;
	Rhs m_rhs;

public:
	couple_product_right(C &&parent, Rhs &&rhs) :
		m_couple(std::forward<C>(parent)),
		m_rhs(std::forward<Rhs>(rhs))
	{
	}

	auto get_first(void) const -> decltype(m_couple.get_first() * m_rhs)
	{
		return m_couple.get_first() * m_rhs;
	}

	auto get_second(void) const -> decltype(m_couple.get_second() * m_rhs)
	{
		return m_couple.get_second() * m_rhs;
	}
};


template <class Lhs, class C>
class couple_product_left
{
private:
	Lhs m_lhs;
	C m_couple;

public:
	couple_product_left(Lhs &&lhs, C &&parent) :
		m_lhs(std::forward<Lhs>(lhs)),
		m_couple(std::forward<C>(parent))
	{
	}

	auto get_first(void) const -> decltype(m_lhs * m_couple.get_first())
	{
		return m_lhs * m_couple.get_first();
	}

	auto get_second(void) const -> decltype(m_lhs * m_couple.get_second())
	{
		return m_lhs * m_couple.get_second();
	}
};


template <class A>
struct is_couple_impl : std::false_type {};

template <class T, class S>
struct is_couple_impl<couple<T, S> > : std::true_type {};

template <class A, class B>
struct is_couple_impl<couple_product_left<A, B> > : std::true_type {};

template <class A, class B>
struct is_couple_impl<couple_product_right<A, B> > : std::true_type {};

template <class A, class B>
struct is_couple_impl<couple_product<A, B> > : std::true_type {};

template <class A>
struct is_couple : is_couple_impl<typename std::decay<A>::type> {};


template <class A, class B>
typename enable_if<
	couple_product_right<A, B>,
	is_couple<A>::value && !is_couple<B>::value
>::type
	operator*(A &&a, B &&b)
{
	return couple_product_right<A, B>(std::forward<A>(a), std::forward<B>(b));
}


template <class A, class B>
typename enable_if<
	couple_product_left<A, B>,
	!is_couple<A>::value && is_couple<B>::value
>::type
	operator*(A &&a, B &&b)
{
	return couple_product_left<A, B>(std::forward<A>(a), std::forward<B>(b));
}


template <class A, class B>
typename enable_if<
	couple_product<A, B>,
	is_couple<A>::value && is_couple<B>::value
>::type
	operator*(A &&a, B &&b)
{
	return couple_product<A, B>(std::forward<A>(a), std::forward<B>(b));
}


template <class C>
typename enable_if<std::ostream &, is_couple<C>::value>::type
	operator<<(std::ostream &os, C const &c)
{
	return os << c.get_first() << ", " << c.get_second();
}


template <class T, class S>
	couple<T, S> create_couple(T &&left, S &&right)
{
	return couple<T, S>(std::forward<T>(left), std::forward<S>(right));
}

#endif // NEWCOUPLE_HPP_INCLUDED

