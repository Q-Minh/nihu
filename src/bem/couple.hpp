/**
 * \file couple.hpp
 * \brief declaration of class ::couple and couple expressions
 * \author Peter Fiala fiala@hit.bme.hu, Peter Rucz rucz@hit.bme.hu
 */

#ifndef COUPLE_HPP_INCLUDED
#define COUPLE_HPP_INCLUDED

#include "../util/product_type.hpp"

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
	typedef typename couple_traits<Derived>::first_value_t first_value_t;
	typedef typename couple_traits<Derived>::second_value_t second_value_t;
	/** \brief type of the first member */
	typedef typename couple_traits<Derived>::first_expr_t first_expr_t;
	/** \brief type of the second member */
	typedef typename couple_traits<Derived>::second_expr_t second_expr_t;

	/**
	 * \brief interface function to return first member
	 * \return const reference to first member
	 */
	first_expr_t first(void) const
	{
		return derived().first();
	}

	/**
	 * \brief interface function to return second member
	 * \return const reference to secon member
	 */
	second_expr_t second(void) const
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
	typedef First first_value_t;
	typedef Second second_value_t;
	typedef First const &first_expr_t;
	typedef Second const &second_expr_t;
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
	typedef typename base_t::first_value_t first_value_t;
	typedef typename base_t::second_value_t second_value_t;
	typedef typename base_t::first_expr_t first_expr_t;
	typedef typename base_t::second_expr_t second_expr_t;

	couple(first_value_t const &f = first_value_t(), second_value_t const &s = second_value_t())
		: m_first(f), m_second(s)
	{
	}

	template <class OtherDerived>
	couple(couple_base<OtherDerived> const &other)
		: m_first(other.first()), m_second(other.second())
	{
	}

	template <class OtherDerived>
	couple &operator=(couple_base<OtherDerived> const &other)
	{
		m_first = other.first();
		m_second = other.second();
		return *this;
	}

	first_expr_t first(void) const
	{
		return m_first;
	}

	first_value_t &first(void)
	{
		return m_first;
	}

	second_expr_t const &second(void) const
	{
		return m_second;
	}

	second_value_t &second(void)
	{
		return m_second;
	}

	template <class OtherDerived>
	couple &operator+=(couple_base<OtherDerived> const &other)
	{
		m_first += other.first();
		m_second += other.second();
		return *this;
	}

protected:
	first_value_t m_first;	/**< \brief the first stored object */
	second_value_t m_second;	/**< \brief the second stored object */
};


template <class LDerived, class R>
struct couple_traits<couple_product_right<LDerived, R> >
{
	typedef typename product_type<typename LDerived::first_value_t, R>::type first_value_t;
	typedef typename product_type<typename LDerived::second_value_t, R>::type second_value_t;
	typedef first_value_t first_expr_t;
	typedef second_value_t second_expr_t;
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
	typedef typename base_t::first_value_t first_value_t;
	typedef typename base_t::second_value_t second_value_t;
	typedef typename base_t::first_expr_t first_expr_t;
	typedef typename base_t::second_expr_t second_expr_t;

	/**
	 * \brief constructor from two term references
	 * \param left reference to left hand side term
	 * \param right reference to left hand side term
	 */
	couple_product_right(LDerived const &left, Right const &right)
		: m_left(left), m_right(right)
	{
	}

	/**
	 * \brief return first object of the product
	 * \return first object
	 */
	first_expr_t first(void) const
	{
		return m_left.first() * m_right;
	}

	/**
	 * \brief return second object of the product
	 * \return second object
	 */
	second_expr_t second(void) const
	{
		return m_left.second() * m_right;
	}

protected:
	LDerived const &m_left;	/**< \brief left hand side term */
	Right const &m_right;	/**< \brief right hand side term */
};


template <class L, class RDerived>
struct couple_traits<couple_product_left<L, RDerived> >
{
	typedef typename product_type<L, typename RDerived::first_value_t>::type first_value_t;
	typedef typename product_type<L, typename RDerived::second_value_t>::type second_value_t;
	typedef first_value_t first_expr_t;
	typedef second_value_t second_expr_t;
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
	typedef typename base_t::first_expr_t first_expr_t;
	typedef typename base_t::second_expr_t second_expr_t;

	/**
	 * \brief constructor from two term references
	 * \param left reference to left hand side term
	 * \param right reference to left hand side term
	 */
	couple_product_left(Left const &left, RDerived const &right)
		: m_left(left), m_right(right)
	{
	}

	/**
	 * \brief return first object of the product
	 * \return first object
	 */
	first_expr_t first(void) const
	{
		return m_left * m_right.first();
	}

	/**
	 * \brief return second object of the product
	 * \return second object
	 */
	second_expr_t second(void) const
	{
		return m_left * m_right.second();
	}

protected:
	Left const &m_left;			/**< \brief left hand side term */
	RDerived const &m_right;	/**< \brief right hand side term */
};

#endif //  COUPLE_HPP_INCLUDED

