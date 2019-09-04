#ifndef P2P_INDEXED_HPP_INCLUDED
#define P2P_INDEXED_HPP_INCLUDED

#include "fmm_operator.hpp"
#include "util/matrix_traits.hpp"
#include <type_traits>

namespace NiHu
{
namespace fmm
{

template <class Operator, class TestIt, class TrialIt>
class p2p_indexed
	: public fmm_operator<p2p_tag>
{
public:
	typedef typename std::decay<Operator>::type operator_t;
	typedef TestIt test_iterator_t;
	typedef TrialIt trial_iterator_t;

	typedef size_t test_input_t;
	typedef size_t trial_input_t;
	typedef typename operator_t::result_t result_t;

	static size_t const num_dof_per_src = num_cols<result_t>::value;
	static size_t const num_dof_per_rec = num_rows<result_t>::value;

	p2p_indexed(Operator &&op, TestIt test_begin, TestIt test_end,
		TrialIt trial_begin, TrialIt trial_end)
		: m_op(std::forward<Operator>(op))
		, m_test_begin(test_begin)
		, m_test_end(test_end)
		, m_trial_begin(trial_begin)
		, m_trial_end(trial_end)
	{
	}

	result_t operator()(size_t test_idx, size_t trial_idx) const
	{
		return m_op(m_test_begin[test_idx], m_trial_begin[trial_idx]);
	}

	size_t num_rec() const
	{
		return m_test_end - m_test_begin;
	}

	size_t num_src() const
	{
		return m_trial_end - m_trial_begin;
	}

private:
	Operator m_op;
	test_iterator_t m_test_begin;
	test_iterator_t m_test_end;
	trial_iterator_t m_trial_begin;
	trial_iterator_t m_trial_end;
};

template <class Operator, class TestIt, class TrialIt>
auto create_p2p_indexed(Operator &&op, TestIt test_begin, TestIt test_end,
	TrialIt trial_begin, TrialIt trial_end)
{
	return p2p_indexed<Operator, TestIt, TrialIt>(
		std::forward<Operator>(op),
		test_begin, test_end,
		trial_begin, trial_end);
}

} // end of namespace fmm
} // namespace NiHu

#endif // P2P_INDEXED_HPP_INCLUDED
