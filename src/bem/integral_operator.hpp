/** \file integral_operator.hpp
* \brief declaration of class ::integral_operator
* \ingroup intop
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
*/

#ifndef INTEGRAL_OPERATOR_HPP_INCLUDED
#define INTEGRAL_OPERATOR_HPP_INCLUDED

#include "single_integral.hpp"
#include "double_integral.hpp"

template <class Scalar, class Derived>
class scaled_integral_operator;

template <class Derived>
struct integral_operator_traits;

template <class Derived>
class integral_operator_base
{
public:
	Derived const &derived() const
	{
		return static_cast<Derived const &>(*this);
	}

	Derived &derived()
	{
		return static_cast<Derived &>(*this);
	}

	typedef integral_operator_traits<Derived> traits_t;

	template <class Test, class Trial>
	struct wr_result_type
	{
		typedef typename traits_t::template wr_result_type<Test, Trial>::type type;
	};

	template <class Test, class Trial>
	typename wr_result_type<Test, Trial>::type
		eval_on_fields(Test const &test, Trial const &trial) const
	{
		return derived().eval_on_fields(test, trial);
	}

/*
	template <class Scalar>
	scaled_integral_operator<Scalar, Derived> operator*(Scalar const &scalar) const
	{
		return ;
	}
*/
};


template <class Scalar, class Derived>
scaled_integral_operator<Scalar, Derived>
	operator*(Scalar const &scalar, integral_operator_base<Derived> const &rhs)
{
	return scaled_integral_operator<Scalar, Derived>(scalar, rhs.derived());
}


template <class Scalar, class Derived>
struct integral_operator_traits<scaled_integral_operator<Scalar, Derived> >
{
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

template <class Scalar, class Derived>
class scaled_integral_operator : public integral_operator_base<scaled_integral_operator<Scalar, Derived> >
{
public:
	typedef integral_operator_base<scaled_integral_operator<Scalar, Derived> > base_t;


	scaled_integral_operator(Scalar const &scalar, integral_operator_base<Derived> const &parent) :
		m_scalar(scalar), m_parent(parent.derived())
	{
	}

	template <class Test, class Trial>
	typename base_t::template wr_result_type<Test, Trial>::type
		eval_on_fields(Test const &test, Trial const &trial) const
	{
		return m_scalar * m_parent.eval_on_fields(test, trial);
	}

private:
	Scalar const m_scalar;
	Derived const m_parent;
};

class identity_integral_operator;

template <>
struct integral_operator_traits<identity_integral_operator>
{
	template <class Test, class Trial>
	struct wr_result_type
	{
		typedef typename single_integral<Test, Trial>::result_t type;
	};
};

class identity_integral_operator : public integral_operator_base<identity_integral_operator>
{
public:
	typedef integral_operator_base<identity_integral_operator> base_t;

	template <class Test, class Trial>
	typename base_t::template wr_result_type<Test, Trial>::type
		eval_on_fields(Test const &test, Trial const &trial) const
	{
		return single_integral<Test, Trial>::eval(test, trial);
	}
};



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



template <class Kernel>
class integral_operator : public integral_operator_base<integral_operator<Kernel> >
{
public:
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

