/**
 * \file couple.hpp
 * \brief declaration of class ::couple and couple expressions
 * \author Peter Fiala fiala@hit.bme.hu, Peter Rucz rucz@hit.bme.hu
 */

#ifndef COUPLE_HPP_INCLUDED
#define COUPLE_HPP_INCLUDED


// forward declaration
template <class Derived>
struct couple_traits;

// forward declaration
template <class LDerived, class Right>
class couple_product_right;

// forward declaration
template <class Left, class RDerived>
class couple_product_left;

/**
 * \brief base class of all couple expressions
 * \tparam Derived CRTP derived class
 */
template <class Derived>
class couple_base
{
private:
	/**
	 * \brief helper function to return reference to CRTP derived class
	 * \return constant reference to derived class
	 */
	Derived const &derived(void) const 
	{
		return static_cast<Derived const&>(*this);
	}

	/**
	 * \brief helper function to return reference to CRTP derived class
	 * \return reference to derived class
	 */
	Derived &derived(void)
	{
		return static_cast<Derived &>(*this);
	}

public:
	/** \brief type of the first member */
	typedef typename couple_traits<Derived>::first_t first_t;
	/** \brief type of the second member */
	typedef typename couple_traits<Derived>::second_t second_t;

	/**
	 * \brief interface function to return first member
	 * \return const reference to first member
	 */
	first_t const &first(void) const
	{
		return derived().first();
	}

	/**
	 * \brief interface function to return second member
	 * \return const reference to secon member
	 */
	second_t const &second(void) const
	{
		return derived().second();
	}

	/**
	 * \brief multiply a couple expression from the right with an arbitrary type
	 * \param rhs the right hand side value
	 * \return a product proxy
	 */
	template <class Right>
	couple_product_right<Derived, Right> operator*(Right const &rhs) const
	{
		return couple_product_right<Derived, Right>(derived(), rhs);
	}

	template <class Left, class RDerived>
	friend couple_product_left<Left, RDerived> operator*(Left const &, couple_base<RDerived> const &);
};

/**
 * \brief multiply a couple expression from the left with an arbitrary type
 * \param l the left hand side factor
 * \param r the right hand side factor
 * \return a product proxy
 */
template <class Left, class RDerived>
couple_product_left<Left, RDerived> operator*(Left const &l, couple_base<RDerived> const &r)
{
	return couple_product_left<Left, RDerived>(l, r.derived());
}


// forward declaration
template <class First, class Second>
class couple;

template <class First, class Second>
struct couple_traits<couple<First, Second> >
{
	typedef First first_t;
	typedef Second second_t;
};

/**
 * \brief class to store two objects of different type
 * \tparam First the type of the first object
 * \tparam Second the type of the second object
 * \details the class is used to conveniently work with kernels that provide
 * two different values
 */
template <class First, class Second = First>
class couple : public couple_base<couple<First, Second> >
{
public:
	typedef couple_base<couple<First, Second> > base_t;
	typedef typename base_t::first_t first_t;
	typedef typename base_t::second_t second_t;
	
	couple(first_t const &f = first_t(), Second const &s = Second()) : m_first(f), m_second(s)
	{
	}
	
	first_t const &first(void) const
	{
		return m_first;
	}
	
	first_t &first(void)
	{
		return m_first;
	}
	
	second_t const &second(void) const
	{
		return m_second;
	}
	
	second_t &second(void)
	{
		return m_second;
	}

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
	first_t m_first;
	second_t m_second;
};


template <class LDerived, class R>
struct couple_traits<couple_product_right<LDerived, R> >
{
	typedef typename LDerived::first_t first_t;
	typedef typename LDerived::second_t second_t;
};

/**
 * \brief class to represent a product proxy of a couple and an arbitrary type
 * \tparam LDerived the couple type
 * \tparam Right the type of the right hand side
 */
template <class LDerived, class Right>
class couple_product_right : public couple_base<couple_product_right<LDerived, Right> >
{
public:
	typedef couple_base<couple_product_right<LDerived, Right> > base_t;
	typedef typename base_t::first_t first_t;
	typedef typename base_t::second_t second_t;

	couple_product_right(LDerived const &left, Right const &right) : m_left(left), m_right(right) {}

	first_t first(void) const
	{
		return m_left.first() * m_right;
	}

	second_t second(void) const
	{
		return m_left.second() * m_right;
	}

protected:
	LDerived const &m_left;
	Right const &m_right;
};


template <class L, class RDerived>
struct couple_traits<couple_product_left<L, RDerived> >
{
	typedef typename RDerived::first_t first_t;
	typedef typename RDerived::second_t second_t;
};

/**
 * \brief class to represent a product proxy of an arbitrary type and a couple
 * \tparam Left the type of the left hand side
 * \tparam RDerived the right hand side couple type
 */
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


#endif //  COUPLE_HPP_INCLUDED
