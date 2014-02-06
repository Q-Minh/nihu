#ifndef CONDITIONAL_PRECOMPUTE_HPP_INCLUDED
#define CONDITIONAL_PRECOMPUTE_HPP_INCLUDED

#include "core/global_definitions.hpp"
#include <type_traits>

<<<<<<< HEAD
/** \brief definition of complexity tags of matrix functions */
namespace matrix_function_complexity
{
	/** \brief indicates that the matrix is a zero expression */
	struct zero { typedef zero type; };
	/** \brief indicates that the matrix is constant and can stored */
	struct constant { typedef constant type; };
	/** \brief indicates that the matrix should be computed on the fly */
=======
namespace matrix_function_complexity
{
	struct zero { typedef zero type; };
	struct constant { typedef constant type; };
>>>>>>> 2e166a6b007b6e4703214e5e6df6377e5769a759
	struct general { typedef general type; };
}


/** \brief Conditionally precompute and store objects
<<<<<<< HEAD
 * \tparam Complexity the result matrix's complexity
 * \tparam Func the class that computes the matrix
 * \tparam Args Arguments passed to the Func class
=======
 * \tparam Complexity true if the object does not need to be precomputed
 * \tparam Func the class that computes object instances
 * \tparam Arguments passed to the Func class
>>>>>>> 2e166a6b007b6e4703214e5e6df6377e5769a759
 */
template <class Complexity, class Func, class ...Args>
class conditional_precompute
{
public:
	/** \brief the return type of the function class */
	typedef decltype(Func::eval(Args()...)) functor_ret_type;
	/** \brief the return type of the static function eval */
	typedef functor_ret_type return_type;

	/** \brief return object computed by the function class */
	static return_type eval(Args const &...args)
	{
		return Func::eval(args...);
	}
};


/** \brief specialisation of ::conditional_precompute for the constant case */
template <class Func, class...Args>
class conditional_precompute<matrix_function_complexity::constant, Func, Args...>
{
public:
	/** \brief the return type of the Func function class */
	typedef decltype(Func::eval(Args()...)) functor_ret_type;
	/** \brief the return type of the static function eval */
	typedef typename std::add_lvalue_reference<
        typename std::add_const<functor_ret_type>::type
    >::type return_type;

	/** \brief return object computed by the function class */
	static return_type eval(Args const &...args)
	{
		static const functor_ret_type m_stored = Func::eval(args...);
		return m_stored;
	}
};


/** \brief specialisation of ::conditional_precompute for the zero case */
template <class Func, class...Args>
class conditional_precompute<matrix_function_complexity::zero, Func, Args...>
{
public:
	/** \brief the return type of the Func function class */
	typedef decltype(Func::eval(Args()...)) functor_ret_type;
	/** \brief the return type of the static function eval */
	typedef decltype(functor_ret_type::Zero()) return_type;

	/** \brief return object computed by the function class */
	static return_type eval(Args const &...args)
	{
		return functor_ret_type::Zero();
	}
};


template <class Complexity, class Func, class ...Args>
class conditional_precompute_instance
{
public:
	/** \brief the return type of the Func function class */
	typedef decltype(Func::eval(Args()...)) functor_ret_type;
	/** \brief the return type of the static function eval */
	typedef functor_ret_type return_type;

	template <class...ConstrArgs>
	conditional_precompute_instance(ConstrArgs const &...)
	{
	}

	/** \brief return object computed by the function class */
	return_type operator()(Args const &...args) const
	{
		return Func::eval(args...);
	}
};


#define ALLOW_PRECOMPUTE

#ifdef ALLOW_PRECOMPUTE

template <class Func, class...Args>
<<<<<<< HEAD
=======
class conditional_precompute<matrix_function_complexity::constant, Func, Args...>
{
public:
	/** \brief the return type of the Func function class */
	typedef decltype(Func::eval(Args()...)) functor_ret_type;
	/** \brief the return type of the static function eval */
	typedef typename std::add_lvalue_reference<
        typename std::add_const<functor_ret_type>::type
    >::type return_type;

	/** \brief return object computed by the function class */
	static return_type eval(Args const &...args)
	{
		static const functor_ret_type m_stored = Func::eval(args...);
		return m_stored;
	}
};


template <class Func, class...Args>
class conditional_precompute<matrix_function_complexity::zero, Func, Args...>
{
public:
	/** \brief the return type of the Func function class */
	typedef decltype(Func::eval(Args()...)) functor_ret_type;
	/** \brief the return type of the static function eval */
	typedef decltype(functor_ret_type::Zero()) return_type;

	/** \brief return object computed by the function class */
	static return_type eval(Args const &...args)
	{
		return functor_ret_type::Zero();
	}
};


template <class Complexity, class Func, class ...Args>
class conditional_precompute_instance
{
public:
	/** \brief the return type of the Func function class */
	typedef decltype(Func::eval(Args()...)) functor_ret_type;
	/** \brief the return type of the static function eval */
	typedef functor_ret_type return_type;

	template <class...ConstrArgs>
	conditional_precompute_instance(ConstrArgs const &...)
	{
	}

	/** \brief return object computed by the function class */
	return_type operator()(Args const &...args) const
	{
		return Func::eval(args...);
	}
};


template <class Func, class...Args>
>>>>>>> 2e166a6b007b6e4703214e5e6df6377e5769a759
class conditional_precompute_instance<matrix_function_complexity::constant, Func, Args...>
{
public:
	/** \brief the return type of the Func function class */
	typedef decltype(Func::eval(Args()...)) functor_ret_type;
	/** \brief the return type of the static function eval */
	typedef typename std::add_lvalue_reference<
        typename std::add_const<functor_ret_type>::type
    >::type return_type;

	conditional_precompute_instance(Args const &...args)
        : m_stored(Func::eval(args...))
	{
	}

	/** \brief return object computed by the function class */
	template <class...EvalArgs>
	return_type operator()(EvalArgs const &...) const
	{
		return m_stored;
	}

private:
    const functor_ret_type m_stored;
<<<<<<< HEAD
};


template <class Func, class...Args>
class conditional_precompute_instance<matrix_function_complexity::zero, Func, Args...>
{
public:
	/** \brief the return type of the Func function class */
	typedef decltype(Func::eval(Args()...)) functor_ret_type;
	/** \brief the return type of the static function eval */
	typedef decltype(functor_ret_type::Zero()) return_type;

	conditional_precompute_instance(Args const &...args)
	{
	}

	/** \brief return object computed by the function class */
	template <class...EvalArgs>
	return_type operator()(EvalArgs const &...) const
	{
		return functor_ret_type::Zero();
	}
=======
>>>>>>> 2e166a6b007b6e4703214e5e6df6377e5769a759
};

#endif

template <class Func, class...Args>
class conditional_precompute_instance<matrix_function_complexity::zero, Func, Args...>
{
public:
	/** \brief the return type of the Func function class */
	typedef decltype(Func::eval(Args()...)) functor_ret_type;
	/** \brief the return type of the static function eval */
	typedef decltype(functor_ret_type::Zero()) return_type;

	conditional_precompute_instance(Args const &...args)
	{
	}

	/** \brief return object computed by the function class */
	template <class...EvalArgs>
	return_type operator()(EvalArgs const &...) const
	{
		return functor_ret_type::Zero();
	}
};

#endif // CONDITIONAL_PRECOMPUTE_HPP_INCLUDED
