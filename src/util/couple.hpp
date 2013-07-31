/**
 * \file couple.hpp
 * \brief declaration of couple expressions
 * \author Peter Fiala fiala@hit.bme.hu, Peter Rucz rucz@hit.bme.hu
 */

#ifndef COUPLE_HPP_INCLUDED
#define COUPLE_HPP_INCLUDED

#include "crtp_base.hpp"
#include "product_type.hpp"
#include <tuple>

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
	 * \return member expression
	 */
	template <size_t idx>
	auto get(void) const
		-> decltype(const_crtp_ptr<Derived, std::integral_constant<int,idx> >()->template get<idx>())
	{
		return derived().template get<idx>();
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
	* \param parent the couple expression whos row is expressed
	* \param idx the row index
	*/
	couple_row(couple &parent, unsigned idx)
		: m_parent(parent), m_idx(idx)
	{
	}

private:
	template <class otherDerived, size_t idx>
	struct increase
	{
		static void eval(
			couple_row const & This,
			couple_base<otherDerived> const &other)
		{
			This.m_parent.template get<idx-1>().row(This.m_idx)
				+= other.template get<idx-1>();
			increase<otherDerived, idx-1>::eval(This, other);
		}
	};

	template <class otherDerived>
	struct increase<otherDerived, 0>
	{
		static void eval(couple_row const &, couple_base<otherDerived> const &)
		{
		}
	};

public:
	/**
	* \brief incremenet a row by an other couple
	* \tparam otherDerived the other couple expression type
	* \param other constant reference to the other couple
	*/
	template <class otherDerived>
	couple_row const &operator += (couple_base<otherDerived> const &other)
	{
		increase<otherDerived, couple::tuple_size>::eval(*this, other);
		return *this;
	}

protected:
	/** \brief reference to the parent couple */
	couple &m_parent;
	/** \brief the row index */
	unsigned const m_idx;
};



template <class...Args>
class couple;

template <class...Args>
struct couple_traits<couple<Args...> >
{
	typedef std::tuple<Args...> tuple_t;
	enum { tuple_size = std::tuple_size<tuple_t>::value };
};


/**
 * \brief class to store a couple
 * \tparam Args the type of the object
 */
template <class...Args>
class couple :
	public couple_base<couple<Args...> >,
	public std::tuple<Args...>
{
public:
	/** \brief self-returning metafunction */
	typedef couple type;
	typedef std::tuple<Args...> base_t;

	typedef couple_traits<type> traits_t;
	enum { tuple_size = traits_t::tuple_size };

 	/** \brief constructor initialising all members
	 * \param [in] args the arguments
	 */
	couple(Args const &...args) :
		base_t(args...)
	{
	}

 	/** \brief constructor initialising all members
	 * \tparam OtherArgs the argument types
	 * \param [in] otherArgs the arguments
	 */
	template <class...OtherArgs>
	couple(OtherArgs&&...otherArgs) :
		base_t(std::forward<OtherArgs>(otherArgs)...)
	{
	}

private:
	template <size_t idx>
	void setZeroImpl(std::integral_constant<size_t, idx>)
	{
		std::get<idx>(*this).setZero();
		setZeroImpl(std::integral_constant<size_t, idx+1>());
	}

	void setZeroImpl(
		std::integral_constant<size_t, tuple_size>)
	{
	}

public:
	/** \brief set both matrices to zero
	 * \todo this function forges ::couple to Eigen matrices. Sick.
	 */
	couple const &setZero(void)
	{
		setZeroImpl(std::integral_constant<size_t, 0>());
		return *this;
	}

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

	template <class C, class OtherDerived, size_t idx>
	struct add_impl
	{
		static void eval(C &c, couple_base<OtherDerived> const &other)
		{
			std::get<idx-1>(c) += other.template get<idx-1>();
			add_impl<C, OtherDerived, idx-1>::eval(c, other);
		}
	};

	template <class C, class OtherDerived>
	struct add_impl<C, OtherDerived, 0>
	{
		static void eval(C &, couple_base<OtherDerived> const &)
		{
		}
	};


public:
 	/** \brief increment with an other couple
	 * \param [in] other the other couple to increment with
	 * \return reference to this
	 */
	template <class OtherDerived>
	couple &operator+=(couple_base<OtherDerived> const &other)
	{
		add_impl<couple, OtherDerived, tuple_size>::eval(*this, other);
		return *this;
	}

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

template <class A, class B>
couple<
	typename std::enable_if<is_eigen<A>::value, A>::type,
	typename std::enable_if<is_eigen<B>::value, B>::type
>
operator,(A &&a, B &&b)
{
	return couple<A, B>(std::forward<A>(a), std::forward<B>(b));
}





template <class LDerived, class Right>
class couple_product_right;


template <class LDerived, class Right>
struct couple_traits<couple_product_right<LDerived, Right> >
{
	template <class Elem>
	struct prod_t
	{
		typedef typename product_type<Elem, Right>::type type;
	};

	template <class C>
	struct tuple2args;

	template <class...Args>
	struct tuple2args<std::tuple<Args...> >
	{
		typedef std::tuple<typename prod_t<Args>::type...> type;
	};

	typedef typename tuple2args<
		typename couple_traits<typename std::decay<LDerived>::type>::tuple_t
	>::type tuple_t;

	enum {
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

	LDerived m_left;	/**< \brief left hand side term */
	Right m_right;		/**< \brief right hand side term */

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




template <class Left, class RDerived>
class couple_product_left;


template <class Left, class RDerived>
struct couple_traits<couple_product_left<Left, RDerived> >
{
	template <class Elem>
	struct prod_t : product_type<Left, Elem> {};

	template <class C>
	struct tuple2args;

	template <class...Args>
	struct tuple2args<std::tuple<Args...> >
	{
		typedef std::tuple<prod_t<Args>...> type;
	};

	typedef typename tuple2args<
		typename couple_traits<typename std::decay<RDerived>::type>::tuple_t
	>::type tuple_t;

	enum {
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
	Left m_left;		/**< \brief left hand side term */
	RDerived m_right;	/**< \brief right hand side term */

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
	 * \brief return elemnt of the couple product
	 * \return element of the couple product
	 */
	template <size_t idx>
	auto get(void) const
		-> decltype( m_left * m_right.template get<idx>() )
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
inline typename std::enable_if<
		is_couple<LDerived>::value,
		couple_product_right<LDerived, Right>
	>::type
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
inline typename std::enable_if<
		is_couple<RDerived>::value,
		couple_product_left<Left, RDerived>
	>::type
	operator*(Left &&lhs, RDerived &&rhs)
{
	return couple_product_left<Left, RDerived>(
		std::forward<Left>(lhs),
		std::forward<RDerived>(rhs));
}

#endif //  COUPLE_HPP_INCLUDED

