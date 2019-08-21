/** \file p2x_integral.hpp
 * \brief implementation of class template p2x_integral
 */

#ifndef P2X_INTEGRAL_HPP_INCLUDED
#define P2X_INTEGRAL_HPP_INCLUDED

#include "core/gaussian_quadrature.hpp"
#include "core/field.hpp"

#include "fmm_operator.hpp"

#include "integral_operator_expression.hpp"
#include "jacobian_computer.hpp"
#include "util/matrix_traits.hpp"
#include "util/type2tag.hpp"

#include <type_traits>

namespace NiHu
{
namespace fmm
{

/** \brief integrate a p2x-operator over a trial field
 * \tparam Operator the operator to integrate
 * \tparam TrialField the trial field type
 */
template <class Operator, class TrialField>
class p2x_integral;


template <class Operator, class TrialField>
struct integral_operator_expression_traits<p2x_integral<Operator, TrialField> >
{
	typedef typename std::decay<Operator>::type operator_t;
	typedef typename operator_t::test_input_t test_input_t;
	typedef TrialField trial_input_t;
	typedef Eigen::Matrix<
		typename scalar<typename operator_t::result_t>::type,
		num_rows<typename operator_t::result_t>::value,
		num_cols<typename operator_t::result_t>::value * TrialField::nset_t::num_nodes
	> result_t;
};


template <class Operator, class TrialField>
class p2x_integral
	: public integral_operator_expression<p2x_integral<Operator, TrialField> >
	, public fmm_operator<typename std::decay<Operator>::type::fmm_tag>
{
public:
	typedef typename std::decay<Operator>::type operator_t;
	typedef TrialField trial_field_t;

	typedef typename trial_field_t::elem_t elem_t;
	typedef typename trial_field_t::nset_t nset_t;

	typedef typename operator_t::test_input_t test_input_t;
	typedef trial_field_t trial_input_t;

	typedef typename elem_t::domain_t domain_t;
	typedef typename domain_t::xi_t xi_t;
	typedef NiHu::gaussian_quadrature<domain_t> quadrature_t;

	typedef typename operator_t::result_t op_result_t;
	static int const op_num_cols = num_cols<op_result_t>::value;
	static int const op_num_rows = num_rows<op_result_t>::value;
	static size_t const result_cols = op_num_cols * nset_t::num_nodes;
	typedef Eigen::Matrix<
		typename scalar<typename operator_t::result_t>::type,
		op_num_rows, result_cols
	> result_t;

	p2x_integral(Operator &&op, size_t order)
		: m_op(std::forward<Operator>(op))
		, m_quadrature(unsigned(order))
	{
	}

	size_t rows(test_input_t const &to) const
	{
		return m_op.rows(to);
	}

	size_t cols(trial_input_t const &) const
	{
		return result_cols;
	}

	result_t operator()(test_input_t const &to, trial_field_t const &field) const
	{
		result_t mat = result_t::Zero(rows(to), cols(field));
		elem_t const &elem = field.get_elem();
		for (auto const &q : m_quadrature)
		{
			xi_t const &xi = q.get_xi();
			auto N = nset_t::template eval_shape<0>(xi);

			double w = q.get_w();
			typename operator_t::trial_input_t tri(elem, xi);
			double jac = jacobian_computer<elem_t>::eval(elem, xi);
			auto K = m_op(to, tri);
			for (Eigen::Index i = 0; i < nset_t::num_nodes; ++i)
				mat.block(0, i*op_num_cols, K.rows(), op_num_cols) += 
				K 
				* (N(i, 0) * Eigen::Matrix<double, op_num_cols, op_num_cols>::Identity())
				*(w * jac);
		}
		return mat;
	}

private:
	Operator m_op;
	quadrature_t m_quadrature;
};

template <class Operator, class TrialFieldTag>
p2x_integral<Operator, typename tag2type<TrialFieldTag>::type>
create_p2x_integral(Operator &&op, size_t order, TrialFieldTag)
{
	typedef typename tag2type<TrialFieldTag>::type TrialField;
	return p2x_integral<Operator, TrialField>(std::forward<Operator>(op), order);
}

} // end of namespace fmm
} // namespace NiHu

#endif // P2X_INTEGRAL_HPP_INCLUDED
