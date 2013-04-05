#ifndef COUPLE_HPP_INCLUDED
#define COUPLE_HPP_INCLUDED

template <class Derived>
struct couple_traits;

template <class LDerived, class RDerived>
class couple_product_right;

template <class LDerived, class Right>
class couple_product_right;

template <class Left, class RDerived>
class couple_product_left;

template <class Derived>
class couple_base
{
private:
	Derived const &derived(void) const 
	{
		return static_cast<Derived const&>(*this);
	}

	Derived &derived(void)
	{
		return static_cast<Derived &>(*this);
	}

protected:

public:
	typedef typename couple_traits<Derived>::first_t first_t;
	typedef typename couple_traits<Derived>::second_t second_t;

	first_t const &first(void) const
	{
		return derived().first();
	}

	second_t const &second(void) const
	{
		return derived().second();
	}

	template <class Right>
	couple_product_right<Derived, Right> operator*(Right const &rhs) const
	{
		return couple_product_right<Derived, Right>(derived(), rhs);
	}

	template <class Left, class RDerived>
	friend couple_product_left<Left, RDerived> operator*(Left const &, couple_base<RDerived> const &);
};

template <class Left, class RDerived>
couple_product_left<Left, RDerived> operator*(Left const &l, couple_base<RDerived> const &r)
{
	return couple_product_left<Left, RDerived>(l, r.derived());
}



template <class First, class Second>
class couple;

template <class First, class Second>
struct couple_traits<couple<First, Second> >
{
	typedef First first_t;
	typedef Second second_t;
};

template <class First, class Second = First>
class couple : public couple_base<couple<First, Second> >
{
public:
	couple(First const &f = First(), Second const &s = Second()) : m_first(f), m_second(s) {}
	First const &first(void) const { return m_first; }
	First &first(void) { return m_first; }
	Second const &second(void) const { return m_second; }
	Second &second(void) { return m_second; }

	template <class OtherDerived>
	couple &operator=(couple_base<OtherDerived> const &other)
	{
		m_first = other.first();
		m_second = other.second();
		return *this;
	}

	template <class OtherDerived>
	couple &operator+=(couple_base<OtherDerived> const &other)
	{
		m_first += other.first();
		m_second += other.second();
		return *this;
	}

protected:
	First m_first;
	Second m_second;
};



template <class LDerived, class R>
struct couple_traits<couple_product_right<LDerived, R> >
{
	typedef typename LDerived::first_t first_t;
	typedef typename LDerived::second_t second_t;
};

template <class Derived, class Right>
class couple_product_right : public couple_base<couple_product_right<Derived, Right> >
{
public:
	typedef couple_base<couple_product_right<Derived, Right> > base_t;
	typedef typename base_t::first_t first_t;
	typedef typename base_t::second_t second_t;

	couple_product_right(Derived const &left, Right const &right) : m_left(left), m_right(right) {}

	first_t first(void) const
	{
		return m_left.first() * m_right;
	}

	second_t second(void) const
	{
		return m_left.second() * m_right;
	}

protected:
	Derived const &m_left;
	Right const &m_right;
};







template <class L, class RDerived>
struct couple_traits<couple_product_left<L, RDerived> >
{
	typedef typename RDerived::first_t first_t;
	typedef typename RDerived::second_t second_t;
};

template <class Left, class RDerived>
class couple_product_left : public couple_base<couple_product_left<Left, RDerived> >
{
public:
	typedef couple_base<couple_product_left<Left, RDerived> > base_t;
	typedef typename base_t::first_t first_t;
	typedef typename base_t::second_t second_t;

	couple_product_left(Left const &left, RDerived const &right) : m_left(left), m_right(right) {}

	first_t first(void) const
	{
		return m_left * m_right.first();
	}

	second_t second(void) const
	{
		return m_left * m_right.second();
	}

protected:
	Left const &m_left;
	RDerived const &m_right;
};





/*


template <class LDerived, class RDerived>
struct couple_traits<couple_product_right<LDerived, couple_base<RDerived> > >
{
	typedef typename LDerived::first_t first_t;
	typedef typename LDerived::second_t second_t;
};

template <class LDerived, class RDerived>
class couple_product_right<LDerived, couple_base<RDerived> > : public couple_base<couple_product_right<LDerived, RDerived > >
{
public:
	typedef couple_base<couple_product_right<LDerived, RDerived> > base_t;
	typedef typename base_t::first_t first_t;
	typedef typename base_t::second_t second_t;

	couple_product_right(LDerived const &left, RDerived const &right) : m_left(left), m_right(right) {}

	first_t first(void) const
	{
		return m_left.first() * m_right.first();
	}

	second_t second(void) const
	{
		return m_left.second() * m_right.second();
	}

protected:
	LDerived const &m_left;
	RDerived const &m_right;
};

*/

#endif //  COUPLE_HPP_INCLUDED
