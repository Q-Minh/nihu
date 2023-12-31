// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2019  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2019  Peter Rucz <rucz@hit.bme.hu>
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

/**\file conditional_precompute.hpp
 * \brief Implementation of conditionally precomputed and stored quantities
 */

#ifndef CONDITIONAL_PRECOMPUTE_HPP_INCLUDED
#define CONDITIONAL_PRECOMPUTE_HPP_INCLUDED

#include "nihu/core/global_definitions.hpp"
#include <type_traits>

namespace NiHu
{

/** \brief definition of complexity tags of functions returning matrices */
namespace matrix_function_complexity
{
	/** \brief the returned matrix is a zero expression */
	struct zero { typedef zero type; };
	/** \brief the returned matrix is constant and can be stored */
	struct constant { typedef constant type; };
	/** \brief the returned matrix should be computed on the fly */
	struct general { typedef general type; };
}


/** \brief Conditionally precompute and store objects
 * \tparam Complexity the result matrix's complexity
 * \tparam Func the class whose static method computes the matrix
 * \tparam Args Arguments passed to the eval method of Func
 */
template <class Complexity, class Func, class ...Args>
class conditional_precompute
{
public:
	/** \brief self-returning */
	typedef conditional_precompute type;
	/** \brief the return type of the function class */
	typedef decltype( Func::eval(std::declval<Args>()...) ) functor_ret_type;
	/** \brief the return type of the static function eval */
	typedef functor_ret_type return_type;

	/** \brief return object computed by the function class */
	static return_type eval(Args const &...args)
	{
		return Func::eval(args...);
	}
};


/** \brief specialisation of NiHu::conditional_precompute for the matrix_function_complexity::constant case */
template <class Func, class...Args>
class conditional_precompute<matrix_function_complexity::constant, Func, Args...>
{
public:
	/** \brief self-returning */
	typedef conditional_precompute type;
	/** \brief the return type of the Func function class */
	typedef decltype( Func::eval(std::declval<Args>()...) ) functor_ret_type;
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


/** \brief specialisation of NiHu::conditional_precompute for the matrix_function_complexity::zero case */
template <class Func, class...Args>
class conditional_precompute<matrix_function_complexity::zero, Func, Args...>
{
public:
	/** \brief self-returning */
	typedef conditional_precompute type;
	/** \brief the return type of the Func function class */
	typedef decltype( Func::eval(std::declval<Args>()...) ) functor_ret_type;
	/** \brief the return type of the static function eval */
	typedef decltype(functor_ret_type::Zero()) return_type;

	/** \brief return object computed by the function class */
	static constexpr return_type eval(Args const &...)
	{
		return functor_ret_type::Zero();
	}
};


template <class Complexity, class Func, class ...Args>
class conditional_precompute_instance
{
public:
	/** \brief self-returning */
	typedef conditional_precompute_instance type;
	/** \brief the return type of the Func function class */
	typedef decltype( Func::eval(std::declval<Args>()...) ) functor_ret_type;
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
class conditional_precompute_instance<matrix_function_complexity::constant, Func, Args...>
{
public:
	/** \brief self-returning */
	typedef conditional_precompute_instance type;
	/** \brief the return type of the Func function class */
	typedef decltype( Func::eval(std::declval<Args>()...) ) functor_ret_type;
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
};


template <class Func, class...Args>
class conditional_precompute_instance<matrix_function_complexity::zero, Func, Args...>
{
public:
	/** \brief self-returning */
	typedef conditional_precompute_instance type;
	/** \brief the return type of the Func function class */
	typedef decltype( Func::eval(std::declval<Args>()...) ) functor_ret_type;
	/** \brief the return type of the static function eval */
	typedef decltype(functor_ret_type::Zero()) return_type;

	conditional_precompute_instance(Args const &...args)
	{
	}

	/** \brief return object computed by the function class */
	template <class...EvalArgs>
	return_type constexpr operator()(EvalArgs const &...) const
	{
		return functor_ret_type::Zero();
	}
};

} // end of namespace NiHu

#endif // CONDITIONAL_PRECOMPUTE_HPP_INCLUDED

