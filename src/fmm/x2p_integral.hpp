/** \file x2p_integral.hpp
 * \brief implementation of class template x2p_integral
 */

#ifndef X2P_INTEGRAL_HPP_INCLUDED
#define X2P_INTEGRAL_HPP_INCLUDED

#include "core/gaussian_quadrature.hpp"
#include "core/field.hpp"
#include "fmm_operator.hpp"

#include "integral_operator_expression.hpp"
#include "jacobian_computer.hpp"
#include "util/matrix_traits.hpp"
#include "util/type2tag.hpp"

#include <type_traits> // std::enable_if

namespace NiHu
{
namespace fmm
{

/** \brief integrate an x2p-operator over a test field
 * \tparam Operator the operator to integrate
 * \tparam TestField the test field type
 */
template <class Operator, class TestField>
class x2p_integral;

template <class Operator, class TestField>
struct integral_operator_expression_traits<x2p_integral<Operator, TestField> >
{
	typedef typename std::decay<Operator>::type operator_t;
	typedef TestField test_input_t;
	typedef typename operator_t::trial_input_t trial_input_t;
	typedef Eigen::Matrix<
		typename scalar<typename operator_t::result_t>::type,
		num_rows<typename operator_t::result_t>::value *TestField::nset_t::num_nodes,
		num_cols<typename operator_t::result_t>::value
	> result_t;
};


/** \brief specialisation for integrals with general test fields */
template <class Operator, class TestField>
class x2p_integral
	: public integral_operator_expression<x2p_integral<Operator, TestField> >
	, public fmm_operator<typename std::decay<Operator>::type::fmm_tag>
{
public:
	typedef integral_operator_expression<x2p_integral<Operator, TestField> > base_t;

	typedef typename std::decay<Operator>::type operator_t;
	typedef TestField test_field_t;

	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::result_t result_t;

	typedef typename test_field_t::elem_t elem_t;
	typedef typename test_field_t::nset_t nset_t;

	typedef typename elem_t::domain_t domain_t;
	typedef typename domain_t::xi_t xi_t;
	typedef NiHu::gaussian_quadrature<domain_t> quadrature_t;

	static size_t const result_rows =
		num_rows<typename operator_t::result_t>::value * nset_t::num_nodes;

	x2p_integral(Operator &&op, size_t order)
		: m_op(std::forward<Operator>(op))
		, m_quadrature(unsigned(order))
	{
	}

	size_t rows(test_input_t const &) const
	{
		return result_rows;
	}

	size_t cols(trial_input_t const &from) const
	{
		return m_op.cols(from);
	}

	/** \brief evaluate the operator between two inputs
	 * \param [in] to the result field
	 * \param [in] from the source cluster or field
	 * The operator evaluates th stored kernel for each integration point over the
	 * field, weights by the weigting function and sums the results
	 */
	result_t operator()(test_input_t const &to, trial_input_t const &from) const
	{
		result_t mat = result_t::Zero(result_rows, cols(from));
		elem_t const &elem = to.get_elem();
		// traverse quadrature elements
		for (auto const &q : m_quadrature)
		{
			xi_t const &xi = q.get_xi();
			double w = q.get_w();
			typename operator_t::test_input_t tsi(elem, xi);
			double jac = jacobian_computer<elem_t>::eval(elem, xi);
			mat += nset_t::template eval_shape<0>(xi) * (w * jac)
				* m_op(tsi, from);
		}
		return mat;
	}

private:
	/** \brief the operator */
	Operator m_op;
	/** \brief the quadrature */
	quadrature_t m_quadrature;
};



/** \brief specialisation for integrals with dirac test fields */
template <class Operator, class TestField>
class x2p_integral<Operator, NiHu::dirac_field<TestField> >
	: public integral_operator_expression<x2p_integral<Operator, NiHu::dirac_field<TestField> > >
	, public fmm_operator<typename std::decay<Operator>::type::fmm_tag>
{
public:
	typedef typename std::decay<Operator>::type operator_t;
	typedef TestField test_field_t;

	typedef integral_operator_expression<x2p_integral<Operator, NiHu::dirac_field<TestField> > > base_t;

	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::result_t result_t;

	typedef typename test_field_t::elem_t elem_t;
	typedef typename test_field_t::nset_t nset_t;

	typedef typename elem_t::domain_t domain_t;
	typedef typename domain_t::xi_t xi_t;
	typedef NiHu::gaussian_quadrature<domain_t> quadrature_t;

	static size_t const result_rows =
		num_rows<typename operator_t::result_t>::value * nset_t::num_nodes;

	x2p_integral(Operator &&op, size_t order = 0)
		: m_op(std::forward<Operator>(op))
	{
	}

	size_t rows(test_input_t const &) const
	{
		return result_rows;
	}

	size_t cols(trial_input_t const &from) const
	{
		return m_op.cols(from);
	}

	result_t operator()(test_input_t const &to, trial_input_t const &from) const
	{
		result_t mat = result_t::Zero(result_rows, cols(from));
		auto const M = num_rows<typename operator_t::result_t>::value;
		for (size_t N = 0; N < nset_t::num_nodes; ++N)
		{
			// get the node in intrinsic coordinates
			xi_t const &xi = nset_t::corner_at(N);
			// form the test input
			typename operator_t::test_input_t tsi(to.get_elem(), xi);
			// compute the operator and place in result's block
			mat.block(N * M, 0, M, cols(from)) = m_op(tsi, from);
		}

		return mat;
	}

private:
	/** \todo should the operator be stored by value? */
	Operator m_op;
};

template <class Operator, class FieldTag>
x2p_integral<Operator, typename tag2type<FieldTag>::type>
create_x2p_integral(Operator &&op, size_t order, FieldTag)
{
	typedef typename tag2type<FieldTag>::type test_field_t;
	return x2p_integral<Operator, test_field_t>(std::forward<Operator>(op), order);
}

} // end of namespace fmm
} // namespace NiHu

#endif // X2P_INTEGRAL_HPP_INCLUDED
