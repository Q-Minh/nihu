/** \file integral_operator.hpp
 * \brief definition of linear arithmetics of integral operators
 */

#ifndef FMM_INTEGRAL_OPERATOR_HPP_INCLUDED
#define FMM_INTEGRAL_OPERATOR_HPP_INCLUDED

#include <utility>
#include <type_traits>

namespace NiHu
{
namespace fmm
{

/** \brief the base class of every integral operator
 * \tparam Derived the CRTP derived class
 */
template <class Derived>
class integral_operator_expression;

/** \brief the traits structure of an integral operator
 * \tparam Derived the CRTP derived class
 */
template <class Derived>
struct integral_operator_expression_traits;

/** \brief metafunction to determine if C is an integral operator
 * \tparam C the input
 */
template <class C>
struct is_integral_operator_expression
	: public std::is_base_of<
	integral_operator_expression<typename std::decay<C>::type>,
	typename std::decay<C>::type> {};


/** \brief the sum of two integral operators
 * \tparam LhsDerived the left hand side type
 * \tparam RhsDerived the right hand side type
 */
template <class LhsDerived, class RhsDerived>
class integral_operator_sum;

/** \brief the difference of two integral operators
 * \tparam LhsDerived the left hand side type
 * \tparam RhsDerived the right hand side type
 */
template <class LhsDerived, class RhsDerived>
class integral_operator_diff;

/** \brief constant times an integral operator
 * \tparam LhsDerived the left hand side type
 * \tparam Scalar the scalar type
 */
template <class LhsDerived, class Scalar>
class integral_operator_scaled;

/** \brief source-concatenation of two integral operators
 * \tparam LhsDerived the left hand side type
 * \tparam LhsDerived the right hand side type
 */
template <class LhsDerived, class RhsDerived>
class integral_operator_src_concatenated;

template <class Derived>
class integral_operator_expression
{
public:
	typedef Derived derived_t;
	typedef integral_operator_expression_traits<derived_t> traits_t;

	typedef typename traits_t::trial_input_t trial_input_t;
	typedef typename traits_t::test_input_t test_input_t;
	typedef typename traits_t::result_t result_t;

	derived_t const &derived() const
	{
		return *(static_cast<derived_t const *>(this));
	}

	derived_t &derived()
	{
		return *(static_cast<derived_t *>(this));
	}

	size_t rows(test_input_t const &ti) const
	{
		return derived().rows(ti);
	}

	size_t cols(trial_input_t const &ti) const
	{
		return derived().cols(ti);
	}

	result_t operator()(test_input_t const &tsi, trial_input_t const &tri) const
	{
		return derived()(tsi, tri);
	}
};

template <class Lhs, class Rhs, typename std::enable_if<
	is_integral_operator_expression<Lhs>::value &
	is_integral_operator_expression<Rhs>::value, int>::type = 0>
integral_operator_sum<Lhs, Rhs> operator+(Lhs &&lhs, Rhs &&rhs)
{
	return integral_operator_sum<Lhs, Rhs>(
		std::forward<Lhs>(lhs),	std::forward<Rhs>(rhs));
}


template <class Lhs, class Rhs, typename std::enable_if<
	is_integral_operator_expression<Lhs>::value &
	is_integral_operator_expression<Rhs>::value, int>::type = 0>
integral_operator_diff<Lhs, Rhs> operator-(Lhs &&lhs, Rhs &&rhs)
{
	return integral_operator_diff<Lhs, Rhs>(
		std::forward<Lhs>(lhs), std::forward<Rhs>(rhs));
}


template <class Lhs, class Scalar,
	typename std::enable_if<is_integral_operator_expression<Lhs>::value, int>::type = 0>
integral_operator_scaled<Lhs, Scalar> operator*(Lhs &&lhs, Scalar &&c)
{
	return integral_operator_scaled<Lhs, Scalar>(
		std::forward<Lhs>(lhs), std::forward<Scalar>(c));
}

template <class Scalar, class Rhs,
	typename std::enable_if<is_integral_operator_expression<Rhs>::value, int>::type = 0>
integral_operator_scaled<Rhs, Scalar>
operator *(Scalar &&c, Rhs &&rhs)
{
	return integral_operator_scaled<Rhs, Scalar>(
		std::forward<Rhs>(rhs), std::forward<Scalar>(c));
}																					  


/// \brief integral operator_expression_traits for the sum of two integral operators
/// \tparam Lhs the left hand side operator
/// \tparam Rhs the right hand side operator
template <class Lhs, class Rhs>
struct integral_operator_expression_traits<integral_operator_sum<Lhs, Rhs> >
{
	typedef typename std::decay<Lhs>::type lhs_derived_t;
	typedef typename integral_operator_expression_traits<lhs_derived_t>::test_input_t test_input_t;
	typedef typename integral_operator_expression_traits<lhs_derived_t>::trial_input_t trial_input_t;
	typedef typename integral_operator_expression_traits<lhs_derived_t>::result_t result_t;
};


template <class Lhs, class Rhs>
struct integral_operator_expression_traits<integral_operator_diff<Lhs, Rhs> >
{
	typedef typename std::decay<Lhs>::type lhs_derived_t;
	typedef typename integral_operator_expression_traits<lhs_derived_t>::test_input_t test_input_t;
	typedef typename integral_operator_expression_traits<lhs_derived_t>::trial_input_t trial_input_t;
	typedef typename integral_operator_expression_traits<lhs_derived_t>::result_t result_t;
};


template <class Lhs, class Rhs>
class integral_operator_sum
	: public integral_operator_expression<integral_operator_sum<Lhs, Rhs> >
{
public:
	typedef integral_operator_expression<integral_operator_sum<Lhs, Rhs> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;

	integral_operator_sum(Lhs &&lhs, Rhs &&rhs)
		: m_lhs(std::forward<Lhs>(lhs))
		, m_rhs(std::forward<Rhs>(rhs))
	{
	}

	size_t rows(test_input_t const &ti) const
	{
		return m_lhs.rows(ti);
	}

	size_t cols(trial_input_t const &ti) const
	{
		return m_lhs.cols(ti);
	}

	result_t operator()(test_input_t const &tsi, trial_input_t const &tri) const
	{
		return m_lhs(tsi, tri) + m_rhs(tsi, tri);
	}

private:
	Lhs m_lhs;
	Rhs m_rhs;
};


template <class Lhs, class Rhs>
class integral_operator_diff
	: public integral_operator_expression<integral_operator_diff<Lhs, Rhs> >
{
public:
	typedef integral_operator_expression<integral_operator_diff<Lhs, Rhs> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;

	integral_operator_diff(Lhs &&lhs, Rhs &&rhs)
		: m_lhs(std::forward<Lhs>(lhs))
		, m_rhs(std::forward<Rhs>(rhs))
	{
	}

	size_t rows(test_input_t const &ti) const
	{
		return m_lhs.rows(ti);
	}

	size_t cols(trial_input_t const &ti) const
	{
		return m_lhs.cols(ti);
	}

	result_t operator()(test_input_t const &tsi, trial_input_t const &tri) const
	{
		return m_lhs(tsi, tri) - m_rhs(tsi, tri);
	}

private:
	Lhs m_lhs;
	Rhs m_rhs;
};


/// \brief specialisation of integral_operator_expression_traits for the scaled integral operator
/// \tparam Lhs the original integral operator
/// \tparam Scalar the scalar type
template <class Lhs, class Scalar>
struct integral_operator_expression_traits<integral_operator_scaled<Lhs, Scalar> >
{
	typedef typename std::decay<Lhs>::type derived_t;
	typedef typename integral_operator_expression_traits<derived_t>::test_input_t test_input_t;
	typedef typename integral_operator_expression_traits<derived_t>::trial_input_t trial_input_t;
	/// @todo should be replaced by product_type (result_t Scalar)
	typedef typename integral_operator_expression_traits<derived_t>::result_t result_t;
};


template <class Lhs, class Scalar>
class integral_operator_scaled
	: public integral_operator_expression<integral_operator_scaled<Lhs, Scalar> >
{
public:
	typedef integral_operator_expression<integral_operator_scaled<Lhs, Scalar> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;

	integral_operator_scaled(Lhs &&lhs, Scalar &&c)
		: m_lhs(std::forward<Lhs>(lhs))
		, m_c(std::forward<Scalar>(c))
	{
	}

	size_t rows(test_input_t const &ti) const
	{
		return m_lhs.rows(ti);
	}

	size_t cols(trial_input_t const &ti) const
	{
		return m_lhs.cols(ti);
	}

	result_t operator()(test_input_t const &tsi, trial_input_t const &tri) const
	{
		return m_lhs(tsi, tri) * m_c;
	}

private:
	Lhs m_lhs;
	Scalar m_c;
};


template <class Lhs, class Rhs>
struct integral_operator_expression_traits<integral_operator_src_concatenated<Lhs, Rhs> >
{
	typedef typename std::decay<Lhs>::type lhs_derived_t;
	typedef typename std::decay<Rhs>::type rhs_derived_t;
	static int const lhs_cols = num_cols<typename lhs_derived_t::result_t>::value;
	static int const lhs_rows = num_rows<typename lhs_derived_t::result_t>::value;
	static int const rhs_cols = num_cols<typename rhs_derived_t::result_t>::value;
	typedef typename scalar<typename lhs_derived_t::result_t>::type scalar_t;
	typedef Eigen::Matrix<scalar_t, lhs_rows, lhs_cols + rhs_cols> result_t;
	typedef typename integral_operator_expression_traits<lhs_derived_t>::test_input_t test_input_t;
	typedef typename integral_operator_expression_traits<lhs_derived_t>::trial_input_t trial_input_t;
};


template <class Lhs, class Rhs>
class integral_operator_src_concatenated
	: public integral_operator_expression<integral_operator_src_concatenated<Lhs, Rhs> >
{
public:
	typedef integral_operator_expression<integral_operator_src_concatenated<Lhs, Rhs> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;

	integral_operator_src_concatenated(Lhs &&lhs, Rhs &&rhs)
		: m_lhs(std::forward<Lhs>(lhs))
		, m_rhs(std::forward<Rhs>(rhs))
	{
	}

	size_t rows(test_input_t const &ti) const
	{
		return m_lhs.rows(ti);
	}

	size_t cols(trial_input_t const &ti) const
	{
		return m_lhs.cols(ti) + m_rhs.cols(ti);
	}

	result_t operator()(test_input_t const &tsi, trial_input_t const &tri) const
	{
		result_t res(rows(tsi), cols(tri));
		res.leftCols(m_lhs.cols(tri)) = m_lhs(tsi, tri);
		res.rightCols(m_rhs.cols(tri)) = m_rhs(tsi, tri);
		return res;
	}

private:
	Lhs m_lhs;
	Rhs m_rhs;
};

template <class Lhs, class Rhs>
integral_operator_src_concatenated<Lhs, Rhs>
src_concatenate(Lhs &&lhs, Rhs &&rhs)
{
	return integral_operator_src_concatenated<Lhs, Rhs>(std::forward<Lhs>(lhs), std::forward<Rhs>(rhs));
}


} // end of namespace fmm
} // namespace NiHu

#endif // FMM_integral_operator_OPERATOR_HPP_INCLUDED
