#ifndef FMBEM_HPP_INCLUDED
#define FMBEM_HPP_INCLUDED

#include "util/type2tag.hpp"
#include "p2p_integral.hpp"
#include "p2x_integral.hpp"
#include "x2p_integral.hpp"

namespace NiHu
{
namespace fmm
{

template <class Fmm, class TestField, class TrialField>
class fmbem
{
public:
	typedef typename std::decay<Fmm>::type fmm_t;
	typedef TestField test_field_t;
	typedef TrialField trial_field_t;

	fmbem(Fmm &&fmm, size_t order, bool singular_check_needed)
		: m_fmm(std::forward<Fmm>(fmm))
		, m_far_field_quadrature_order(order)
		, m_singular_check_needed(singular_check_needed)
	{
	}


	template <int Ny>
	p2x_integral<typename fmm_t::template p2m<Ny>, trial_field_t>
		create_p2m() const
	{
		return create_p2x_integral(m_fmm.template create_p2m<Ny>(),
			m_far_field_quadrature_order,
			typename type2tag<trial_field_t>::type());
	}

	template <int Ny>
	p2x_integral<typename fmm_t::template p2l<Ny>, trial_field_t>
		create_p2l() const
	{
		return create_p2x_integral(m_fmm.template create_p2l<Ny>(),
			m_far_field_quadrature_order,
			typename type2tag<trial_field_t>::type());
	}

	template <int Nx>
	x2p_integral<typename fmm_t::template m2p<Nx>, test_field_t>
		create_m2p() const
	{
		return create_x2p_integral(m_fmm.template create_m2p<Nx>(),
			m_far_field_quadrature_order,
			typename type2tag<test_field_t>::type());
	}

	template <int Nx>
	x2p_integral<typename fmm_t::template l2p<Nx>, test_field_t>
		create_l2p() const
	{
		return create_x2p_integral(m_fmm.template create_l2p<Nx>(),
			m_far_field_quadrature_order,
			typename type2tag<test_field_t>::type());
	}

	template <int Nx, int Ny>
	p2p_integral<
		typename fmm_t::template p2p_type<Nx, Ny>::type,
		test_field_t,
		trial_field_t
	>
		create_p2p() const
	{
		return create_p2p_integral(
			m_fmm.template create_p2p<Nx, Ny>(),
			typename type2tag<test_field_t>::type(),
			typename type2tag<trial_field_t>::type(),
			m_singular_check_needed);
	}

private:
	Fmm m_fmm;
	size_t m_far_field_quadrature_order;
	bool m_singular_check_needed;
};

template <class Fmm, class TestTag, class TrialTag>
fmbem<Fmm, typename tag2type<TestTag>::type, typename tag2type<TrialTag>::type>
create_fmbem(Fmm &&fmm, TestTag, TrialTag, size_t far_field_quadrature_order, bool sing_needed)
{
	return fmbem<
		Fmm,
		typename tag2type<TestTag>::type,
		typename tag2type<TrialTag>::type
	>(std::forward<Fmm>(fmm), far_field_quadrature_order, sing_needed);
}

} // namespace fmm
} // namespace NiHu

#endif // FMBEM_HPP_INCLUDED