#ifndef FMM_INTEGRATED_HPP_INCLUDED
#define FMM_INTEGRATED_HPP_INCLUDED

#include "p2x_integral.hpp"
#include "x2p_integral.hpp"
#include "p2p_integral.hpp"


namespace NiHu
{
namespace fmm
{

template <class Op, class FmmTag = typename std::decay<Op>::type::fmm_tag>
class integrated;

template <class Op>
class integrated<Op, m2p_tag> {
public:
	template <class TestTag, class TrialTag>
	static auto eval(Op &&op, TestTag test_tag, TrialTag trial_tag, size_t quadrature_order, bool sing_check)
	{
		return create_x2p_integral(
				std::forward<Op>(op), quadrature_order, test_tag);
	}
};


template <class Op>
class integrated<Op, l2p_tag> {
public:
	template <class TestTag, class TrialTag>
	static auto eval(Op &&op, TestTag test_tag, TrialTag trial_tag, size_t quadrature_order, bool sing_check)
	{
		return create_x2p_integral(
			std::forward<Op>(op), quadrature_order, test_tag);
	}
};


template <class Op>
class integrated<Op, p2m_tag> {
public:
	template <class TestTag, class TrialTag>
	static auto eval(Op &&op, TestTag test_tag, TrialTag trial_tag, size_t quadrature_order, bool sing_check)
	{
		return create_p2x_integral(
			std::forward<Op>(op), quadrature_order, trial_tag);
	}
};


template <class Op>
class integrated<Op, p2l_tag> {
public:
	template <class TestTag, class TrialTag>
	static auto eval(Op &&op, TestTag test_tag, TrialTag trial_tag, size_t quadrature_order, bool sing_check)
	{
		return create_p2x_integral(
			std::forward<Op>(op), quadrature_order, trial_tag);
	}
};


template <class Op>
class integrated<Op, p2p_tag> {
public:
	template <class TestTag, class TrialTag>
	static auto eval(Op &&op, TestTag test_tag, TrialTag trial_tag, size_t quadrature_order, bool sing_check)
	{
		return create_p2p_integral(std::forward<Op>(op), test_tag, trial_tag, sing_check);
	}
};



template <class Op, class TestTag, class TrialTag>
auto create_integrated(Op &&op, TestTag test_tag, TrialTag trial_tag, size_t quadrature_order, bool sing_check)
{
	return integrated<Op>::eval(std::forward<Op>(op), test_tag, trial_tag, quadrature_order, sing_check);
}


template <class TestTag, class TrialTag>
struct integrated_functor
{
	integrated_functor(size_t quadrature_order, bool sing_check)
		: m_quadrature_order(quadrature_order)
		, m_sing_check(sing_check)
	{
	}

	template <class Op>
	auto operator()(Op &&op) const
	{
		return create_integrated(std::forward<Op>(op), TestTag(), TrialTag(),
			m_quadrature_order,
			m_sing_check);
	}

	size_t m_quadrature_order;
	bool m_sing_check;
};

template <class TestTag, class TrialTag>
auto create_integrated_functor(
	TestTag test_tag, TrialTag trial_tag,
	size_t quadrature_order, bool sing_check)
{
	return integrated_functor<TestTag, TrialTag>(
		quadrature_order, sing_check);
}

}	// end of namespace fmm
}	// end of namespace NiHu

#endif /* FMM_INDEXED_HPP_INCLUDED */