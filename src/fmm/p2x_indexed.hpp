/// \file p2x_indexed.hpp
/// \brief implementation of class fmm::p2x_indexed
#ifndef P2X_INDEXED_HPP_INCLUDED
#define P2X_INDEXED_HPP_INCLUDED

#include "fmm_operator.hpp"
#include <type_traits>

namespace NiHu
{
namespace fmm
{

/// \brief class performing indexing of a P2P operator
/// \tparam Operator the original P2P operator
/// \tparam It the indexing iterator type
template <class Operator, class It>
class p2x_indexed
	: public fmm_operator<typename std::decay<Operator>::type::fmm_tag>

{
public:
	typedef typename std::decay<Operator>::type operator_t;
	typedef It iterator_t;

	typedef typename operator_t::test_input_t test_input_t;
	typedef size_t trial_input_t;
	typedef typename operator_t::result_t result_t;

	p2x_indexed(Operator &&op, It begin, It end)
		: m_op(std::forward<Operator>(op))
		, m_begin(begin)
		, m_end(end)
	{
	}

	result_t operator()(test_input_t const &ti, size_t idx) const
	{
		return m_op(ti, typename operator_t::trial_input_t(m_begin[idx]));
	}

	template <class TSI>
	result_t operator()(TSI const &ti, size_t idx) const
	{
		return m_op(test_input_t(ti), typename operator_t::trial_input_t(m_begin[idx]));
	}

	size_t num_src() const
	{
		return m_end - m_begin;
	}

private:
	Operator m_op;
	iterator_t m_begin;
	iterator_t m_end;
};

template <class Operator, class It>
p2x_indexed<Operator, It>
create_p2x_indexed(Operator &&op, It begin, It end)
{
	return p2x_indexed<Operator, It>(std::forward<Operator>(op), begin, end);
}

} // end of namespace fmm
} // namespace NiHu

#endif // P2X_INDEXED_HPP_INCLUDED
