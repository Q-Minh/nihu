// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

/**
 * \file couple.hpp
 * \brief Couple expressions
 */

#ifndef COUPLE_HPP_INCLUDED
#define COUPLE_HPP_INCLUDED

#include "../core/global_definitions.hpp"

#include "crtp_base.hpp"
#include "product_type.hpp"
#include <tuple>

/** \brief traits class of couple expressions */
template <class C>
struct couple_traits;

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
	 * \brief interface function to return member expression
	 * \tparam idx the expression index
	 * \return member expression
	 */
	template <size_t idx>
	auto get(void) const
		-> decltype(const_crtp_ptr<Derived, std::integral_constant<int,idx> >()->template get<idx>())
	{
		return derived().template get<idx>();
	}

private:
	/** \brief specialisation of ::eval_to_tuple for the 1-argument tuple case */
	typename couple_traits<Derived>::tuple_t
	eval_to_tuple_impl(std::integral_constant<unsigned, 1>) const
	{
		return typename couple_traits<Derived>::tuple_t(get<0>());
	}

	/** \brief specialisation of ::eval_to_tuple for the 2-argument tuple case */
	typename couple_traits<Derived>::tuple_t
	eval_to_tuple_impl(std::integral_constant<unsigned, 2>) const
	{
		return typename couple_traits<Derived>::tuple_t(get<0>(), get<1>());
	}

	/** \brief specialisation of ::eval_to_tuple for the 3-argument tuple case */
	typename couple_traits<Derived>::tuple_t
	eval_to_tuple_impl(std::integral_constant<unsigned, 3>) const
	{
		return typename couple_traits<Derived>::tuple_t(get<0>(), get<1>(), get<2>());
	}

	/** \brief specialisation of ::eval_to_tuple for the 4-argument tuple case */
	typename couple_traits<Derived>::tuple_t
	eval_to_tuple_impl(std::integral_constant<unsigned, 4>) const
	{
		return typename couple_traits<Derived>::tuple_t(get<0>(), get<1>(), get<2>(), get<3>());
	}

public:
	/** \brief evaluate the expression into a tuple value class */
	typename couple_traits<Derived>::tuple_t
	eval_to_tuple(void) const
	{
		/** \brief parameter dispatching */
		return eval_to_tuple_impl(
			std::integral_constant<unsigned, couple_traits<Derived>::tuple_size>());
	}
};


/** \brief a row of a couple of matrices
 * \tparam couple the parent couple type
 * \todo This is a minimal implementation. couple_row should be derived from couple_base
 */
template <class couple>
class couple_row
{
public:
	/** \brief constructor
	* \param [in] parent the couple expression whose row is expressed
	* \param [in] row_idx the row index
	*/
	couple_row(couple &parent, unsigned row_idx)
		: m_parent(parent), m_idx(row_idx)
	{
	}

	/** \brief recursive implementation of operator += */
	template <class otherDerived, size_t idx>
	struct increase
	{
		/** \brief evaluate couple row
		 * \param [in] This couple expression
		 * \param [in] other couple expression
		 */
		static void eval(
			couple_row const & This,
			couple_base<otherDerived> const &other)
		{
			This.m_parent.template get<idx-1>().row(This.m_idx)
				+= other.template get<idx-1>();
			increase<otherDerived, idx-1>::eval(This, other);
		}
	};

	/** \brief terminating case of increase */
	template <class otherDerived>
	struct increase<otherDerived, 0>
	{
		/** \brief evaluate empty couple */
		static void eval(couple_row const &, couple_base<otherDerived> const &)
		{
		}
	};

public:
	/**
	* \brief increment a row by a couple expression
	* \tparam otherDerived the other couple expression type
	* \param other constant reference to the other couple
	*/
	template <class otherDerived>
	couple_row const &operator += (couple_base<otherDerived> const &other)
	{
		increase<otherDerived, couple_traits<couple>::tuple_size>::eval(*this, other);
		return *this;
	}

protected:
	/** \brief reference to the parent couple */
	couple &m_parent;
	/** \brief the row index */
	unsigned const m_idx;
};


// forward declaration
template <class...Args>
class couple;

/** \brief traits of a couple */
template <class...Args>
struct couple_traits<couple<Args...> >
{
	/** \brief the underlying tuple class */
	typedef std::tuple<Args...> tuple_t;
	/** \brief the tuple size */
	enum { tuple_size = std::tuple_size<tuple_t>::value };
};


/**
 * \brief class to store a couple
 * \tparam Args the type vector of the object
 */
template <class...Args>
class couple :
	public couple_base<couple<Args...> >,
	public std::tuple<Args...>
{
public:
	/** \brief self-returning metafunction */
	typedef couple type;
	/** \brief the tuple base */
	typedef std::tuple<Args...> base_t;

	/** \brief the traits class */
	typedef couple_traits<type> traits_t;

	/** \brief default constructor */
	couple() {}

 	/** \brief constructor initialising all members
	 * \param [in] args the arguments
	 */
	couple(Args const &...args) :
		base_t(args...)
	{
	}

	/** \brief constructor from a couple expression */
	template <class OtherDerived>
	couple(couple_base<OtherDerived> const &other) :
		base_t(other.eval_to_tuple())
	{
	}

	base_t eval_to_tuple(void) const
	{
		return *this;
	}

private:
	/** \brief recursive implementation of setZero */
	template <size_t idx>
	void setZeroImpl(std::integral_constant<size_t, idx>)
	{
		std::get<idx-1>(*this).setZero();
		setZeroImpl(std::integral_constant<size_t, idx-1>());
	}

	/** \brief terminating case of setZero */
	void setZeroImpl(std::integral_constant<size_t, 0>) {}

public:
	/** \brief set all member matrices zero */
	couple const &setZero(void)
	{
		setZeroImpl(std::integral_constant<size_t, traits_t::tuple_size>());
		return *this;
	}

	/** \brief return a zero couple */
	static couple Zero(void)
	{
		couple res;
		return res.setZero();
	}

 	/** \brief return couple element
	 * \return couple element
	 */
	template <size_t idx>
	typename std::tuple_element<idx, base_t>::type const &
		get(void) const
	{
		return std::get<idx>(*this);
	}

 	/** \brief return reference to couple element
	 * \return reference to couple element
	 */
	template <size_t idx>
	typename std::tuple_element<idx, std::tuple<Args...> >::type &
		get(void)
	{
		return std::get<idx>(*this);
	}

private:
	/** \brief recursive implementation of operator += */
	template <class C, class OtherDerived, size_t idx>
	struct add_impl
	{
		static void eval(C &c, couple_base<OtherDerived> const &other)
		{
			std::get<idx-1>(c) += other.template get<idx-1>();
			add_impl<C, OtherDerived, idx-1>::eval(c, other);
		}
	};

	/** \brief terminating case of operator += */
	template <class C, class OtherDerived>
	struct add_impl<C, OtherDerived, 0>
	{
		static void eval(C &, couple_base<OtherDerived> const &) {}
	};


public:
 	/** \brief increment with an other couple
	 * \param [in] other the other couple to increment with
	 * \return reference to this
	 */
	template <class OtherDerived>
	couple &operator+=(couple_base<OtherDerived> const &other)
	{
		add_impl<couple, OtherDerived, traits_t::tuple_size>::eval(*this, other);
		return *this;
	}

	/** \brief return row proxy */
	couple_row<couple> row(unsigned idx)
	{
		return couple_row<couple>(*this, idx);
	}
};

/** \brief factory function of a couple class */
template <class...Args>
couple<Args...> create_couple(Args &&...args)
{
	return couple<Args...>(std::forward<Args>(args)...);
}


// forward declaration
template <class LDerived, class Right>
class couple_product_right;

/** \brief traits class of a couple multiplied from the right */
template <class LDerived, class Right>
struct couple_traits<couple_product_right<LDerived, Right> >
{
	/** \brief a helper metafunction */
	template <class Elem>
	struct prod_t : product_type<Elem, Right> {};

	/** \brief convert a tuple to a variadic argument pack */
	template <class C>
	struct tuple2args;

	/** \brief convert a tuple to a variadic argument pack */
	template <class...Args>
	struct tuple2args<std::tuple<Args...> >
	{
		typedef std::tuple<typename prod_t<Args>::type...> type;
	};

	/** \brief the tuple type of the product */
	typedef typename tuple2args<
		typename couple_traits<typename std::decay<LDerived>::type>::tuple_t
	>::type tuple_t;

	enum {
		/** \brief the tuple size of the product */
		tuple_size = couple_traits<typename std::decay<LDerived>::type>::tuple_size
	};
};



/**
 * \brief class to represent a product of a couple and an arbitrary type
 * \tparam LDerived the couple type
 * \tparam Right the type of the right hand side
 */
template <class LDerived, class Right>
class couple_product_right :
	public couple_base<couple_product_right<LDerived, Right> >
{
protected:
	/** \brief left hand side term */
	LDerived m_left;
	/** \brief right hand side term */
	Right m_right;

	typedef couple_traits<couple_product_right> traits_t;

public:
	/**
	 * \brief constructor from two term references
	 * \param left reference to left hand side term
	 * \param right reference to left hand side term
	 */
	couple_product_right(LDerived &&left, Right &&right) :
		m_left(std::forward<LDerived>(left)),
		m_right(std::forward<Right>(right))
	{
	}

	/**
	 * \brief return element of the couple product
	 * \return element of the couple product
	 */
	template <size_t idx>
	auto get(void) const
		-> decltype( m_left.template get<idx>() * m_right )
	{
		return m_left.template get<idx>() * m_right;
	}
};


// forward declaration
template <class Left, class RDerived>
class couple_product_left;


/**
 * \brief traits of a product of an arbitrary type and a couple
 * \tparam Left the type of the left hand side
 * \tparam RDerived the right hand side couple type
 */
template <class Left, class RDerived>
struct couple_traits<couple_product_left<Left, RDerived> >
{
	/** \brief helper function */
	template <class Elem>
	struct prod_t : product_type<Left, Elem> {};

	/** \brief general declaration */
	template <class C>
	struct tuple2args;

	/** \brief specialisation for tuples */
	template <class...Args>
	struct tuple2args<std::tuple<Args...> >
	{
		typedef std::tuple<typename prod_t<Args>::type...> type;
	};

	/** \brief the underlying tuple type of the product */
	typedef typename tuple2args<
		typename couple_traits<typename std::decay<RDerived>::type>::tuple_t
	>::type tuple_t;

	enum {
		/** \brief the size of the tuple product */
		tuple_size = couple_traits<typename std::decay<RDerived>::type>::tuple_size
	};
};


/**
 * \brief class to represent a product proxy of an arbitrary type and a couple
 * \tparam Left the type of the left hand side
 * \tparam RDerived the right hand side couple type
 */
template <class Left, class RDerived>
class couple_product_left :
	public couple_base<couple_product_left<Left, RDerived> >
{
protected:
	/** \brief left hand side term */
	Left m_left;
	/** \brief right hand side term */
	RDerived m_right;

public:
	/**
	 * \brief constructor from two term references
	 * \param left reference to left hand side term
	 * \param right reference to left hand side term
	 */
	couple_product_left(Left &&left, RDerived &&right) :
		m_left(std::forward<Left>(left)),
		m_right(std::forward<RDerived>(right))
	{
	}

	/**
	 * \brief return element of the couple product
	 * \return element of the couple product
	 */
	template <size_t idx> auto get(void) const
		-> decltype ( m_left * m_right.template get<idx>() ) /* from GCC4.8 */
		//-> decltype( This->get_left() * This->get_right().template get<idx>() )
	{
		return m_left * m_right.template get<idx>();
	}
};


/** \brief metafunction determining if argument is a couple expression
 * \tparam T the class to investigate
 */
template <class T>
struct is_couple : std::is_base_of<
	couple_base<typename std::decay<T>::type>,
	typename std::decay<T>::type
> {};


/**
 * \brief multiply a couple expression from the right by an arbitrary type
 * \tparam LDerived the left hand side derived type
 * \tparam Right the right hand side type
 * \param [in] lhs the left hand side value
 * \param [in] rhs the right hand side value
 * \return a product proxy
 */
template <class LDerived, class Right>
inline couple_product_right<
	typename std::enable_if<is_couple<LDerived>::value, LDerived>::type,
	Right
>
	operator*(LDerived &&lhs, Right &&rhs)
{
	return couple_product_right<LDerived, Right>(
		std::forward<LDerived>(lhs),
		std::forward<Right>(rhs));
}

/**
 * \brief multiply a couple expression from the left by an arbitrary type
 * \tparam Left the left hand side type
 * \tparam Right the right hand side type
 * \param [in] lhs the left hand side factor
 * \param [in] rhs the right hand side factor
 * \return a product proxy
 */
template <class Left, class RDerived>
inline couple_product_left<
	Left,
	typename std::enable_if<is_couple<RDerived>::value, RDerived>::type
>
	operator*(Left &&lhs, RDerived &&rhs)
{
	return couple_product_left<Left, RDerived>(
		std::forward<Left>(lhs),
		std::forward<RDerived>(rhs));
}

#endif //  COUPLE_HPP_INCLUDED

