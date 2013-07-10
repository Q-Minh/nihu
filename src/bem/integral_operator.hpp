/** \file integral_operator.hpp
* \brief declaration of class ::integral_operator
* \ingroup intop
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
*/

#ifndef INTEGRAL_OPERATOR_HPP_INCLUDED
#define INTEGRAL_OPERATOR_HPP_INCLUDED

#include "single_integral.hpp"
#include "double_integral.hpp"

// forward declaration
template <class Scalar, class Derived>
class scaled_integral_operator;

/* \brief traits class for an integral operator
* \tparam Derived the CRTP derived class
*/
template <class Derived>
struct integral_operator_traits;

/* \brief CRTP base of integral operator expressions
* \tparam Derived the CRTP derived class
*/
template <class Derived>
class integral_operator_base
{
public:
	/** \brief CRTP helper function */
	Derived const &derived() const
	{
		return static_cast<Derived const &>(*this);
	}

	/** \brief CRTP helper function */
	Derived &derived()
	{
		return static_cast<Derived &>(*this);
	}

	/** \brief the traits class of the derived integral operator */
	typedef integral_operator_traits<Derived> traits_t;

	/** \brief metafunction obtained from the traits class */
	template <class Test, class Trial>
	struct wr_result_type
	{
		typedef typename traits_t::template wr_result_type<Test, Trial>::type type;
	};

	/** \brief sub-weighted residual on a test and a trial field
	* \tparam Test the test field type
	* \tpara Trial the trial field type
	* \param [in] test the test field reference
	* \param [in] trial the trial field reference
	* \return the result submatrix
	*/
	template <class Test, class Trial>
	typename wr_result_type<Test, Trial>::type
		eval_on_fields(Test const &test, Trial const &trial) const
	{
		return derived().eval_on_fields(test, trial);
	}
};


/** \brief factory operator to create a scaled integral operator
* \tparam Scalar the scalar type to multply with
* \tparam Derived the integral operator
* \param [in] scalar the scalar to multiply with
* \param [in] rhs the right hand side integral operator
* \return a scaled integral operator proxy object
*/
template <class Scalar, class Derived>
scaled_integral_operator<Scalar, Derived>
	operator*(Scalar const &scalar, integral_operator_base<Derived> const &rhs)
{
	return scaled_integral_operator<Scalar, Derived>(scalar, rhs.derived());
}


/** \brief traits class of class ::scaled_integral_operator */
template <class Scalar, class Derived>
struct integral_operator_traits<scaled_integral_operator<Scalar, Derived> >
{
	/** \brief metafunction returning the result type of a double integral */
	template <class Test, class Trial>
	struct wr_result_type
	{
		typedef typename plain_type<
			typename product_type<
			Scalar,
			typename integral_operator_traits<Derived>::template wr_result_type<Test, Trial>::type
			>::type
		>::type type;
	};
};

/** \brief Proxy class representing an integral operator multiplied by a scalar
* \tparam Scalar the scalar type
* \tparam Derived the integral operator's type
*/
template <class Scalar, class Derived>
class scaled_integral_operator :
	public integral_operator_base<scaled_integral_operator<Scalar, Derived> >
{
public:
	/** \brief the CRTP base class */
	typedef integral_operator_base<scaled_integral_operator<Scalar, Derived> > base_t;


	/** \brief constructor from a scalar and an integral operator instance
	* \param [in] scalar the scalar instance
	* \param parent the integral operator to multiply with the scalar
	*/
	scaled_integral_operator(
		Scalar const &scalar,
		integral_operator_base<Derived> const &parent) :
		m_scalar(scalar), m_parent(parent.derived())
	{
	}

	/** \brief evaluate a scaled integral operator on a test and a trial field
	* \tparam Test the test field's type
	* \tparam Trial the trial field's type
	* \param [in] test the test field
	* \param [in] triel the trial field
	* \return the result matrix of the double integral
	*/
	template <class Test, class Trial>
	typename base_t::template wr_result_type<Test, Trial>::type
		eval_on_fields(Test const &test, Trial const &trial) const
	{
		return m_scalar * m_parent.eval_on_fields(test, trial);
	}

private:
	/** \brief the scalar miltiplier */
	Scalar const m_scalar;
	/** \brief the parent operator */
	Derived const m_parent;
};


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
};

/** \brief The identity integral operator with kernel \f$\delta(x-y)\f$ */
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
	* \param [in] triel the trial field
	* \return the result matrix of the double integral
	*/
	template <class Test, class Trial>
	typename base_t::template wr_result_type<Test, Trial>::type
		eval_on_fields(Test const &test, Trial const &trial) const
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
	typedef Kernel kernel_t;

	/** \brief constructor from kernel reference
	* \param [in] kernel reference to the kernel
	*/
	integral_operator(Kernel const &kernel) :
		m_kernel(kernel)
	{
	}

	/** \brief return kernel reference
	* \return reference to the kernel
	*/
	Kernel &get_kernel(void) const
	{
		return m_kernel;
	}

	/** \brief evaluate an integral operator on a test and a trial field
	* \tparam Test the test field's type
	* \tparam Trial the trial field's type
	* \param [in] test the test field
	* \param [in] triel the trial field
	* \return the result matrix of the double integral
	*/
	template <class Test, class Trial>
	typename base_t::template wr_result_type<Test, Trial>::type
		eval_on_fields(Test const &test, Trial const &trial) const
	{
		return double_integral<kernel_t, Test, Trial>::eval(
			m_kernel,
			test,
			trial);
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
	create_integral_operator(Kernel const &kernel)
{
	return integral_operator<Kernel>(kernel);
}


#endif

