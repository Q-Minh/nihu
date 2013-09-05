/** \file integral_operator.hpp
* \brief declaration of class ::integral_operator
* \ingroup intop
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
*/

#ifndef INTEGRAL_OPERATOR_HPP_INCLUDED
#define INTEGRAL_OPERATOR_HPP_INCLUDED

#include "../util/crtp_base.hpp"
#include "function_space.hpp"
#include "single_integral.hpp"
#include "double_integral.hpp"

template <class Operator, class TrialSpace>
class integral_transform;

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
	template <class Test, class Trial>
	struct wr_result_type : traits_t::template wr_result_type<Test, Trial> {};

	/** \brief sub-weighted residual on a test and a trial field
	* \tparam Test the test field type
	* \tparam Trial the trial field type
	* \param [in] test the test field reference
	* \param [in] trial the trial field reference
	* \return the result submatrix
	*/
	template <class Test, class Trial, class OnSameMesh>
	typename wr_result_type<Test, Trial>::type
		eval_on_fields(
			field_base<Test> const &test,
			field_base<Trial> const &trial,
			OnSameMesh) const
	{
		return derived().derived_eval_on_fields(test, trial, OnSameMesh());
	}

	/** \brief apply the integral operator on a function space and create an ::integral_transform
	* \tparam FuncSpace the trial function space
	* \param [in] funcspace the function space reference
	* \return integral_transform proxy object
	*/
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

template <class IntOp>
struct is_integral_operator : std::is_base_of<
	integral_operator_base<typename std::decay<IntOp>::type>,
	typename std::decay<IntOp>::type
>{};


template <class Scalar, class IntOp>
class scaled_integral_operator;


/** \brief traits class of class ::scaled_integral_operator */
template <class Scalar, class IntOp>
struct integral_operator_traits<scaled_integral_operator<Scalar, IntOp> >
{
	/** \brief metafunction returning the result type of a double integral */
	template <class Test, class Trial>
	struct wr_result_type : plain_type<
		typename product_type<
		Scalar,
		typename integral_operator_traits<
			typename std::decay<IntOp>::type
		>::template wr_result_type<Test, Trial>::type
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
	* \tparam Test the test field's type
	* \tparam Trial the trial field's type
	* \param [in] test the test field
	* \param [in] trial the trial field
	* \return the result matrix of the double integral
	*/
	template <class Test, class Trial, class OnSameMesh>
	typename base_t::template wr_result_type<Test, Trial>::type
		derived_eval_on_fields(
			field_base<Test> const &test,
			field_base<Trial> const &trial,
			OnSameMesh) const
	{
		return m_scalar * m_parent.eval_on_fields(test, trial, OnSameMesh());
	}

private:
	/** \brief the scalar miltiplier */
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
	template <class Test, class Trial>
	struct wr_result_type
	{
		typedef typename single_integral<Test, Trial>::result_t type;
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
	* \tparam Test the test field's type
	* \tparam Trial the trial field's type
	* \param [in] test the test field
	* \param [in] trial the trial field
	* \return the result matrix of the double integral
	*/
	template <class Test, class Trial, class OnSameMesh = std::false_type>
	typename base_t::template wr_result_type<Test, Trial>::type
		derived_eval_on_fields(
			field_base<Test> const &test,
			field_base<Trial> const &trial,
			OnSameMesh) const
	{
		return single_integral<Test, Trial>::eval(test, trial);
	}
};



// forward declaration
template <class Kernel>
class integral_operator;


template <class Kernel>
struct integral_operator_traits<integral_operator<Kernel> >
{
	template <class Test, class Trial>
	struct wr_result_type
	{
		typedef typename double_integral<Kernel, Test, Trial>::result_t type;
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

	/** \brief constructor from kernel reference
	* \param [in] kernel reference to the kernel
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
	* \tparam Test the test field's type
	* \tparam Trial the trial field's type
	* \param [in] test the test field
	* \param [in] trial the trial field
	* \return the result matrix of the double integral
	*/
	template <class Test, class Trial, class OnSameMesh>
	typename base_t::template wr_result_type<Test, Trial>::type
		derived_eval_on_fields(
			field_base<Test> const &test,
			field_base<Trial> const &trial,
			OnSameMesh) const
	{
		return double_integral<kernel_t, Test, Trial>::eval(
			m_kernel, test, trial, OnSameMesh());
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


/** \brief factory function of an integral operator with couple kernels
* \tparam K1 the kernel type
* \tparam Kernels the remaining kernels' type
* \param [in] k1 the first kernel
* \param [in] kernels the remaining kernels
* \return the integral operator object
*/
template <class K1, class...Kernels>
integral_operator<couple_kernel<K1, Kernels...> >
	create_integral_operator(K1 &&k1, Kernels &&...kernels)
{
	return integral_operator<couple_kernel<K1, Kernels...> >(
		create_couple_kernel(std::forward<K1>(k1), std::forward<Kernels>(kernels)...)
	);
}


#endif

