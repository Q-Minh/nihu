#ifndef P2P_HPP_INCLUDED
#define P2P_HPP_INCLUDED

#include "nihu/util/matrix_traits.hpp"

#include "fmm_operator.hpp"

namespace NiHu
{
namespace fmm
{

template <class Kernel>
class p2p
	: public fmm_operator<p2p_tag>
{
public:
	typedef Kernel kernel_t;
	typedef typename kernel_t::test_input_t test_input_t;
	typedef typename kernel_t::trial_input_t trial_input_t;
	typedef typename kernel_t::result_t result_t;

	p2p(kernel_t const &kernel)
		: m_kernel(kernel)
	{
	}

	size_t rows(test_input_t const &) const
	{
		return num_rows<result_t>::value;
	}

	size_t cols(trial_input_t const &) const
	{
		return num_cols<result_t>::value;
	}

	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		return m_kernel(x, y);
	}

	template <class TSI, class TRI>
	result_t operator()(TSI const &tsi, TRI const &tri) const
	{
		return m_kernel(test_input_t(tsi), trial_input_t(tri));
	}

	kernel_t const &get_kernel() const
	{
		return m_kernel;
	}

private:
	kernel_t m_kernel;
};

} // end of namespace fmm
} // namespace NiHu

#endif
