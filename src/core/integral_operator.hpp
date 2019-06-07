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

/// \file integral_operator.hpp
/// \brief declaration of class NiHu::integral_operator
/// \ingroup intop

#ifndef INTEGRAL_OPERATOR_HPP_INCLUDED
#define INTEGRAL_OPERATOR_HPP_INCLUDED

#include "../util/crtp_base.hpp"
#include "function_space.hpp"
#include "single_integral.hpp"
#include "double_integral.hpp"
#include "integral_transform_forward.hpp"

namespace NiHu
{


/** \brief traits class for an integral operator
 * \tparam Derived the CRTP derived class
 */
template <class Derived>
struct integral_operator_traits;

/** \brief CRTP base of integral operator expressions
 * \tparam Derived the CRTP derived class
 */
template <class Derived>
class integral_operator_base
{
public:
	NIHU_CRTP_HELPERS

	/** \brief the traits class of the derived integral operator */
	typedef integral_operator_traits<Derived> traits_t;

	/** \brief metafunction obtained from the traits class */
	template <class TestField, class TrialField>
	struct wr_result_type : traits_t::template wr_result_type<TestField, TrialField> {};

	/** \brief sub-weighted residual on a test and a trial field
	 * \tparam TestField the test field type
	 * \tparam TrialField the trial field type
	 * \param [in] test_field the test field reference
	 * \param [in] trial_field the trial field reference
	 * \return the result submatrix
	 */
	template <class TestField, class TrialField, class OnSameMesh>
	typename wr_result_type<TestField, TrialField>::type
		eval_on_fields(
			field_base<TestField> const &test_field,
			field_base<TrialField> const &trial_field,
			OnSameMesh) const
	{
		return derived().derived_eval_on_fields(test_field, trial_field, OnSameMesh());
	}

	/// \brief apply the integral operator on a function space and create a NiHu::integral_transform
	/// \tparam FuncSpace the trial function space
	/// \param [in] funcspace the function space reference
	/// \return integral_transform proxy object
	template <class FuncSpace>
	integral_transform<
		Derived,
		typename std::enable_if<is_function_space<FuncSpace>::value, FuncSpace>::type
	>
		operator[](FuncSpace &&funcspace)
	{
		return integral_transform<Derived, FuncSpace>(
			derived(),
			std::forward<FuncSpace>(funcspace));
	}
};

/** \brief metafunction returning true if IntOp is an integral operator expression */
template <class IntOp>
struct is_integral_operator : std::is_base_of<
	integral_operator_base<typename std::decay<IntOp>::type>,
	typename std::decay<IntOp>::type
>{};


// forward declaration
template <class Scalar, class IntOp>
class scaled_integral_operator;


/// \brief traits class of class NiHu::scaled_integral_operator
template <class Scalar, class IntOp>
struct integral_operator_traits<scaled_integral_operator<Scalar, IntOp> >
{
	/** \brief metafunction returning the result type of a double integral */
	template <class TestField, class TrialField>
	struct wr_result_type : plain_type<
		typename product_type<
		Scalar,
		typename integral_operator_traits<
			typename std::decay<IntOp>::type
		>::template wr_result_type<TestField, TrialField>::type
		>::type
	> {};

	/** \brief indicates if the operator is to be evaluated only on the same element */
	static bool const is_local = integral_operator_traits<
		typename std::decay<IntOp>::type
	>::is_local;
};


/** \brief Proxy class representing an integral operator multiplied by a scalar
 * \tparam Scalar the scalar type
 * \tparam IntOp the integral operator's type
 */
template <class Scalar, class IntOp>
class scaled_integral_operator :
	public integral_operator_base<scaled_integral_operator<Scalar, IntOp> >
{
public:
	/** \brief the CRTP base class */
	typedef integral_operator_base<scaled_integral_operator<Scalar, IntOp> > base_t;


	/** \brief constructor from a scalar and an integral operator instance
	 * \param [in] scalar the scalar instance
	 * \param parent the integral operator to multiply with the scalar
	 */
	scaled_integral_operator(
		Scalar &&scalar,
		IntOp &&parent) :
		m_scalar(std::forward<Scalar>(scalar)),
		m_parent(std::forward<IntOp>(parent))
	{
	}

	/** \brief evaluate a scaled integral operator on a test and a trial field
	 * \tparam TestField the test field's type
	 * \tparam TrialField the trial field's type
	 * \param [in] test_field the test field
	 * \param [in] trial_field the trial field
	 * \return the result matrix of the double integral
	 */
	template <class TestField, class TrialField, class OnSameMesh>
	typename base_t::template wr_result_type<TestField, TrialField>::type
		derived_eval_on_fields(
			field_base<TestField> const &test_field,
			field_base<TrialField> const &trial_field,
			OnSameMesh) const
	{
		return m_scalar * m_parent.eval_on_fields(test_field, trial_field, OnSameMesh());
	}

private:
	/** \brief the scalar multiplier */
	Scalar m_scalar;
	/** \brief the parent operator */
	IntOp m_parent;
};


/** \brief factory operator to create a scaled integral operator
 * \tparam Scalar the scalar type to multply with
 * \tparam IntOp the integral operator
 * \param [in] scalar the scalar to multiply with
 * \param [in] intop the right hand side integral operator
 * \return a scaled integral operator proxy object
 */
template <class Scalar, class IntOp>
scaled_integral_operator<
	Scalar,
	typename std::enable_if<is_integral_operator<IntOp>::value, IntOp>::type
>
	operator*(Scalar &&scalar, IntOp &&intop)
{
	return scaled_integral_operator<Scalar, IntOp>(
		std::forward<Scalar>(scalar),
		std::forward<IntOp>(intop));
}


// forward declaration
class identity_integral_operator;

/** \brief traits class of the identity integral operator */
template <>
struct integral_operator_traits<identity_integral_operator>
{
	template <class TestField, class TrialField>
	struct wr_result_type
	{
		typedef typename single_integral<TestField, TrialField>::result_t type;
	};

	/** \brief indicates if the operator is to be evaluated only on the same element */
	static bool const is_local = true;
};

/** \brief The identity integral operator \f$ K(x,y) = \delta(x-y) \f$ */
class identity_integral_operator :
	public integral_operator_base<identity_integral_operator>
{
public:
	/** \brief CRTP base */
	typedef integral_operator_base<identity_integral_operator> base_t;

	/** \brief evaluate an identity operator on a test and a trial field
	 * \tparam TestField the test field's type
	 * \tparam TrialField the trial field's type
	 * \param [in] test the test field
	 * \param [in] trial the trial field
	 * \return the result matrix of the double integral
	 */
	template <class TestField, class TrialField, class OnSameMesh = std::false_type>
	typename base_t::template wr_result_type<TestField, TrialField>::type
		derived_eval_on_fields(
			field_base<TestField> const &test_field,
			field_base<TrialField> const &trial_field,
			OnSameMesh) const
	{
		return single_integral<TestField, TrialField>::eval(test_field, trial_field);
	}
};



// forward declaration
template <class Kernel>
class integral_operator;

/** \brief traits of an integral operator
 * \tparam Kernel the kernel type
 */
template <class Kernel>
struct integral_operator_traits<integral_operator<Kernel> >
{
	/** \brief metafunction returning the weighted residual return type */
	template <class TestField, class TrialField>
	struct wr_result_type
	{
		typedef typename double_integral<typename std::decay<Kernel>::type, TestField, TrialField>::result_t type;
	};

	/** \brief indicates if the operator is to be evaluated only on the same element */
	static bool const is_local = false;
};



/** \brief the general integral operator with an arbitrary kernel
 * \tparam the Kernel class
 */
template <class Kernel>
class integral_operator :
	public integral_operator_base<integral_operator<Kernel> >
{
public:
	/** \brief the CRTP base class */
	typedef integral_operator_base<integral_operator<Kernel> > base_t;

	/** \brief template argument as nested type */
	typedef typename std::decay<Kernel>::type kernel_t;

	/** \brief constructor from kernel
	 * \param [in] kernel the kernel
	 */
	integral_operator(Kernel &&kernel) :
		m_kernel(std::forward<Kernel>(kernel))
	{
	}

	/** \brief return kernel (reference)
	 * \return reference to the kernel
	 */
	Kernel get_kernel(void) const
	{
		return m_kernel;
	}

	/** \brief evaluate an integral operator on a test and a trial field
	 * \tparam TestField the test field's type
	 * \tparam TrialField the trial field's type
	 * \param [in] test the test field
	 * \param [in] trial the trial field
	 * \return the result matrix of the double integral
	 */
	template <class TestField, class TrialField, class OnSameMesh>
	typename base_t::template wr_result_type<TestField, TrialField>::type
		derived_eval_on_fields(
			field_base<TestField> const &test_field,
			field_base<TrialField> const &trial_field,
			OnSameMesh) const
	{
		return double_integral<kernel_t, TestField, TrialField>::eval(
			m_kernel, test_field, trial_field, OnSameMesh());
	}

private:
	/** \brief he underlying kernel */
	Kernel m_kernel;
};


/** \brief factory function of an integral operator
 * \tparam Kernel the kernel type
 * \param [in] kernel the kernel
 * \return the integral operator object
 */
template <class Kernel>
integral_operator<Kernel>
	create_integral_operator(Kernel &&kernel)
{
	return integral_operator<Kernel>(std::forward<Kernel>(kernel));
}

} // end of namespace NiHu

#endif

