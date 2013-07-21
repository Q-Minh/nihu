/**
 * \file couple.hpp
 * \brief declaration of class ::couple and couple expressions
 * \author Peter Fiala fiala@hit.bme.hu, Peter Rucz rucz@hit.bme.hu
 */

#ifndef COUPLE_HPP_INCLUDED
#define COUPLE_HPP_INCLUDED

#include "../util/crtp_base.hpp"

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
	/**
	 * \brief interface function to return first member expression
	 * \return first member expression
	 */
	template <class Dummy = void>
	auto first(void) const
		-> decltype(static_cast<typename ignore<Derived, Dummy>::type const*>(nullptr)->first())
	{
		return derived().first();
	}

	/**
	 * \brief interface function to return second member expression
	 * \return second member expression
	 */
	template <class Dummy = void>
	auto second(void) const
		-> decltype(static_cast<typename ignore<Derived, Dummy>::type const*>(nullptr)->second())
	{
		return derived().second();
	}
	

	template <class Dummy = void>
	auto eval(void) const
		-> decltype(static_cast<typename ignore<Derived, Dummy>::type const*>(nullptr)->eval())
	{
		return derived().eval();
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


/**
 * \brief class to store two objects of different type
 * \tparam First the type of the first object
 * \tparam Second the type of the second object
 */
template <class First, class Second = First>
class couple :
	public couple_base<couple<First, Second> >
{
public:
	/** \brief self-returning metafunction */
	typedef couple type;
	typedef First first_t;
	typedef Second second_t;
	
 	/** \brief constructor initialising both members
	 * \param [in] f the first member
	 * \param [in] s the second member
	 */
	couple(First const &f = First(), Second const &s = Second()) :
		m_first(f), m_second(s)
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
	First const &first(void) const
	{
		return m_first;
	}

 	/** \brief return first member
	 * \return reference to first value
	 */
	First &first(void)
	{
		return m_first;
	}

 	/** \brief return second member
	 * \return second member expression
	 */
	Second const &second(void) const
	{
		return m_second;
	}

 	/** \brief return second member
	 * \return reference to second value
	 */
	Second &second(void)
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
	First m_first;	/**< \brief the first stored object */
	Second m_second;	/**< \brief the second stored object */
};


/**
 * \brief class to represent a product expression of a couple and an arbitrary type
 * \tparam LDerived the couple type
 * \tparam Right the type of the right hand side
 */
template <class LDerived, class Right>
class couple_product_right : public couple_base<couple_product_right<LDerived, Right> >
{
protected:
	LDerived const &m_left;	/**< \brief left hand side term */
	Right const &m_right;	/**< \brief right hand side term */

public:
	typedef decltype(m_left.first() * m_right) first_t;
	typedef decltype(m_left.second() * m_right) second_t;
	
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
	auto first(void) const
		-> decltype(m_left.first() * m_right)
	{
		return m_left.first() * m_right;
	}

	/**
	 * \brief return second object of the product
	 * \return second object
	 */
	auto second(void) const
		-> decltype(m_left.second() * m_right)
	{
		return m_left.second() * m_right;
	}
};


/**
 * \brief class to represent a product proxy of an arbitrary type and a couple
 * \tparam Left the type of the left hand side
 * \tparam RDerived the right hand side couple type
 */
template <class Left, class RDerived>
class couple_product_left : public couple_base<couple_product_left<Left, RDerived> >
{
protected:
	Left const &m_left;			/**< \brief left hand side term */
	RDerived const &m_right;	/**< \brief right hand side term */

public:
	typedef couple_base<couple_product_left<Left, RDerived> > base_t;	/**< \brief base type */

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
	auto first(void) const -> decltype(m_left * m_right.first())
	{
		return m_left * m_right.first();
	}

	/**
	 * \brief return second object of the product
	 * \return second object
	 */
	auto second(void) const -> decltype(m_left * m_right.second())
	{
		return m_left * m_right.second();
	}
};


/** \brief factory function of a couple class */
template <class L, class R>
couple<L, R> create_couple(L const &l, R const &r)
{
	return couple<L, R>(l, r);
}

#endif //  COUPLE_HPP_INCLUDED

