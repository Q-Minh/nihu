/**
 * \file couple.hpp
 * \brief declaration of class ::couple and couple expressions
 * \author Peter Fiala fiala@hit.bme.hu, Peter Rucz rucz@hit.bme.hu
 */

#ifndef COUPLE_HPP_INCLUDED
#define COUPLE_HPP_INCLUDED

#include "../util/crtp_base.hpp"
#include "../util/product_type.hpp"

// forward declaration
/** \brief traits class of a couple expression */
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
	NIHU_CRTP_HELPERS

public:
	/** \brief the traits type */
	typedef couple_traits<Derived> traits_t;
	/** \brief value type of the first member */
	typedef typename traits_t::first_value_t first_value_t;
	/** \brief value type of the second member */
	typedef typename traits_t::second_value_t second_value_t;
	/** \brief expression type of the first member */
	typedef typename traits_t::first_expr_t first_expr_t;
	/** \brief expression type of the second member */
	typedef typename traits_t::second_expr_t second_expr_t;

	/**
	 * \brief interface function to return first member expression
	 * \return first member expression
	 */
	first_expr_t first(void) const
	{
		return derived().first();
	}

	/**
	 * \brief interface function to return second member expression
	 * \return second member expression
	 */
	second_expr_t second(void) const
	{
		return derived().second();
	}

	/**
	 * \brief multiply a couple expression from the right with an arbitrary type
	 * \param [in] rhs the right hand side value
	 * \return a product proxy
	 */
	template <class Right>
	couple_product_right<Derived, Right> operator*(Right const &rhs) const
	{
		return couple_product_right<Derived, Right>(derived(), rhs);
	}

	template <class Left, class RDerived>
	friend couple_product_left<Left, RDerived> operator*(
		Left const &,
		couple_base<RDerived> const &);
};

/**
 * \brief multiply a couple expression from the left with an arbitrary type
 * \param [in] lhs the left hand side factor
 * \param [in] rhs the right hand side factor
 * \return a product proxy
 */
template <class Left, class RDerived>
inline couple_product_left<Left, RDerived> operator*(
	Left const &lhs,
	couple_base<RDerived> const &rhs)
{
	return couple_product_left<Left, RDerived>(lhs, rhs.derived());
}


/** \brief a row of a couple of matrices
 * \tparam couple the parent couple type
 * \todo This is a minimal implementation. couple_row should be derived from couple_base
 */
template <class couple>
class couple_row
{
public:
	/** \brief constructor
	* \param parent the couple expression whos row is expressed
	* \param idx the row index
	*/
	couple_row(couple &parent, unsigned idx)
		: m_parent(parent), m_idx(idx)
	{
	}
	
	/**
	* \brief incremenet a row by an other couple
	* \tparam otherDerived the other couple expression type
	* \param other constant reference to the other couple
	*/
	template <class otherDerived>
	couple_row const &operator += (couple_base<otherDerived> const &other)
	{
		m_parent.first().row(m_idx) += other.first();
		m_parent.second().row(m_idx) += other.second();
		return *this;
	}
	
protected:
	/** \brief reference to the parent couple */
	couple &m_parent;
	/** \brief the row index */
	unsigned const m_idx;
};


// forward declaration
template <class First, class Second>
class couple;

/** \brief specialisation of ::couple_traits for a ::couple
 * \tparam First the first type of the couple
 * \tparam Second the second type of the couple
 */
template <class First, class Second>
struct couple_traits<couple<First, Second> >
{
	typedef First first_value_t;	/**< \brief value type of the first member */
	typedef Second second_value_t;	/**< \brief value type of the second member */
	typedef First const &first_expr_t;	/**< \brief expression type of the first member */
	typedef Second const &second_expr_t;	/**< \brief expression type of the second member */
};

/**
 * \brief class to store two objects of different type
 * \tparam First the type of the first object
 * \tparam Second the type of the second object
 */
template <class First, class Second = First>
class couple : public couple_base<couple<First, Second> >
{
public:
	typedef couple_base<couple<First, Second> > base_t;	/**< \brief base type */
	typedef typename base_t::first_value_t first_value_t;	/**< \brief first member's value type */
	typedef typename base_t::second_value_t second_value_t;	/**< \brief second member's value type */
	typedef typename base_t::first_expr_t first_expr_t;	/**< \brief first member's expression type */
	typedef typename base_t::second_expr_t second_expr_t;	/**< \brief second member's expression type */

 	/** \brief constructor initialising both members
	 * \param [in] f the first member
	 * \param [in] s the second member
	 */
	couple(first_value_t const &f = first_value_t(),
		   second_value_t const &s = second_value_t())
		: m_first(f), m_second(s)
	{
	}

	/** \brief set both matrices to zero
	 * \todo this function forges ::couple to Eigen matrices. Sick.
	 */
	couple const &setZero(void)
	{
		m_first.setZero();
		m_second.setZero();
		return *this;
	}

	static couple Zero(void)
	{
		couple res;
		return res.setZero();
	}

 	/** \brief return first member
	 * \return first member expression
	 */
	first_expr_t first(void) const
	{
		return m_first;
	}

 	/** \brief return first member
	 * \return reference to first value
	 */
	first_value_t &first(void)
	{
		return m_first;
	}

 	/** \brief return second member
	 * \return second member expression
	 */
	second_expr_t const &second(void) const
	{
		return m_second;
	}

 	/** \brief return second member
	 * \return reference to second value
	 */
	second_value_t &second(void)
	{
		return m_second;
	}

 	/** \brief increment with an other couple
	 * \param [in] other the other couple to increment with
	 * \return reference to this
	 */
	template <class OtherDerived>
	couple &operator+=(couple_base<OtherDerived> const &other)
	{
		m_first += other.first();
		m_second += other.second();
		return *this;
	}
	
	couple_row<couple> row(unsigned idx)
	{
		return couple_row<couple>(*this, idx);
	}

protected:
	first_value_t m_first;	/**< \brief the first stored object */
	second_value_t m_second;	/**< \brief the second stored object */
};


/**
 * \brief specialisation of ::couple_traits for ::couple_product_right
 * \tparam LDerived the left hand side couple type
 * \tparam R the type of the right hand side
 */
template <class LDerived, class R>
struct couple_traits<couple_product_right<LDerived, R> >
{
	/** \brief the value type of the first member */
	typedef typename product_type<typename LDerived::first_value_t, R>::type first_value_t;
	/** \brief the value type of the second member */
	typedef typename product_type<typename LDerived::second_value_t, R>::type second_value_t;
	/** \brief the expression type of the first member */
	typedef first_value_t first_expr_t;
	/** \brief the expression type of the second member */
	typedef second_value_t second_expr_t;
};

/**
 * \brief class to represent a product expression of a couple and an arbitrary type
 * \tparam LDerived the couple type
 * \tparam Right the type of the right hand side
 */
template <class LDerived, class Right>
class couple_product_right : public couple_base<couple_product_right<LDerived, Right> >
{
public:
	typedef couple_base<couple_product_right<LDerived, Right> > base_t;	/**< \brief base type */
	typedef typename base_t::first_expr_t first_expr_t;	/**< \brief first expression type */
	typedef typename base_t::second_expr_t second_expr_t;	/**< \brief second expression type */

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


/**
 * \brief specialisation of ::couple_traits for ::couple_product_left
 * \tparam L the type of the left hand side
 * \tparam RDerived the right hand side couple type
 */
template <class L, class RDerived>
struct couple_traits<couple_product_left<L, RDerived> >
{
	/** \brief the value type of the first member */
	typedef typename product_type<L, typename RDerived::first_value_t>::type first_value_t;
	/** \brief the value type of the second member */
	typedef typename product_type<L, typename RDerived::second_value_t>::type second_value_t;
	/** \brief the expression type of the first member */
	typedef first_value_t first_expr_t;
	/** \brief the expression type of the second member */
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
	typedef couple_base<couple_product_left<Left, RDerived> > base_t;	/**< \brief base type */
	typedef typename base_t::first_expr_t first_expr_t;	/**< \brief first expression type */
	typedef typename base_t::second_expr_t second_expr_t;	/**< \brief second expression type */

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


/** \brief factory function of a couple class */
template <class L, class R>
couple<L, R> create_couple(L const &l, R const &r)
{
	return couple<L, R>(l, r);
}

#endif //  COUPLE_HPP_INCLUDED

