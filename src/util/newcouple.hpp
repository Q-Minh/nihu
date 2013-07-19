#ifndef NEWCOUPLE_HPP_INCLUDED
#define NEWCOUPLE_HPP_INCLUDED

#include <utility> // pair
#include <iostream> // pair

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


template <class C, class Rhs>
class product_proxy
{
private:
	C m_couple;
	Rhs m_rhs;

public:
	product_proxy(C &&parent, Rhs rhs) :
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


template <class A>
struct is_couple : std::false_type {};

template <class T, class S>
struct is_couple<couple<T, S> > : std::true_type {};

template <class R, bool cond>
struct enable_if { };

template <class R>
struct enable_if<R, true>
{
	typedef R type;
};

template <class A, class B>
typename enable_if<
	product_proxy<A, B>,
	is_couple<typename std::decay<A>::type>::value
>::type
	operator*(A &&a, B &&b)
{
	return product_proxy<A, B>(std::forward<A>(a), std::forward<B>(b));
}


template <class T>
void print_type_info(void)
{
	std::cout << (std::is_const<typename std::remove_reference<T>::type >::value ? "const " : "");
	std::cout << (std::is_same<typename std::decay<T>::type, int>::value ? "int " : "");
	std::cout << (std::is_lvalue_reference<T>::value ? "&" : "");
	std::cout << (std::is_rvalue_reference<T>::value ? "&&" : "");
}

template <class T, class S>
std::ostream &operator<<(std::ostream &os, couple<T, S> const &c)
{
	return os << c.get_first() << ", " << c.get_second();
}

template <class T, class S>
	couple<T, S> create_couple(T &&left, S &&right)
{
	std::cout << "creating couple< ";
	print_type_info<T>();
	std::cout << " , ";
	print_type_info<S>();
	std::cout << ">\n";
	return couple<T, S>(std::forward<T>(left), std::forward<S>(right));
}

#endif // NEWCOUPLE_HPP_INCLUDED

