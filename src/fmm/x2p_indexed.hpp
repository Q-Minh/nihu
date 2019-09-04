#ifndef X2P_INDEXED_HPP_INCLUDED
#define X2P_INDEXED_HPP_INCLUDED

#include "fmm_operator.hpp"
#include <type_traits>

namespace NiHu
{
namespace fmm
{

template <class Operator, class It>
class x2p_indexed
	: public fmm_operator<typename std::decay<Operator>::type::fmm_tag>
{
public:
	typedef typename std::decay<Operator>::type operator_t;
	typedef It iterator_t;

	typedef typename operator_t::trial_input_t trial_input_t;
	typedef size_t test_input_t;
	typedef typename operator_t::result_t result_t;

	x2p_indexed(Operator &&op, It begin, It end)
		: m_op(std::forward<Operator>(op))
		, m_begin(begin)
		, m_end(end)
	{
	}

	result_t operator()(size_t idx, trial_input_t const &ti) const
	{
		return m_op(m_begin[idx], ti);
	}

	size_t num_rec() const
	{
		return m_end - m_begin;
	}

private:
	Operator m_op;
	iterator_t m_begin;
	iterator_t m_end;
};

template <class Operator, class It>
auto create_x2p_indexed(Operator &&op, It begin, It end)
{
	return x2p_indexed<Operator, It>(std::forward<Operator>(op), begin, end);
}

} // end of namespace fmm
} // namespace NiHu

#endif // X2P_INDEXED_HPP_INCLUDED
