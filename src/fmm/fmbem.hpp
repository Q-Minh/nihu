/// \file fmbem.hpp
/// \brief convert an fmm object into an fmbem one

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

/// \brief class performing integration of the particle-related operators
/// \tparam Fmm the original fmm object's type
/// \tparam TestField the test field
/// \tparam TrialField the trial field
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
	
	/// \brief metafunction returning the P2M operator's type
	/// \tparam Ny the order of source differentiation
	template <int Ny>
	struct p2m_type
	{
		typedef p2x_integral<
			typename fmm_t::template p2m_type<Ny>::type,
			trial_field_t
		> type;
	};


	/// \brief metafunction returning the M2P operator's type
	/// \tparam Nx the order of receiver differentiation
	template <int Nx>
	struct m2p_type
	{
		typedef x2p_integral<
			typename fmm_t::template m2p_type<Nx>::type,
			test_field_t
		> type;
	};


	/// \brief metafunction returning the P2L operator's type
	/// \tparam Ny the order of source differentiation
	template <int Ny>
	struct p2l_type
	{
		typedef p2x_integral<
			typename fmm_t::template p2l_type<Ny>::type,
			trial_field_t
		> type;
	};


	/// \brief metafunction returning the L2P operator's type
	/// \tparam Nx the order of receiver differentiation
	template <int Nx>
	struct l2p_type
	{
		typedef x2p_integral<
			typename fmm_t::template l2p_type<Nx>::type,
			test_field_t
		> type;
	};


	/// \brief metafunction returning the P2P operator's type
	/// \tparam Nx the order of receiver differentiation
	/// \tparam Ny the order of source differentiation
	template <int Nx, int Ny>
	struct p2p_type
	{
		typedef p2p_integral<
			typename fmm_t::template p2p_type<Nx, Ny>::type,
			test_field_t,
			trial_field_t
		> type;
	};
	

	/// \brief factory function for the P2M operator
	/// \tparam Ny the order of source differentiation
	template <int Ny>
	typename p2m_type<Ny>::type create_p2m() const
	{
		return create_p2x_integral(m_fmm.template create_p2m<Ny>(),
			m_far_field_quadrature_order,
			typename type2tag<trial_field_t>::type());
	}

	/// \brief factory function for the P2L operator
	/// \tparam Ny the order of source differentiation
	template <int Ny>
	typename p2l_type<Ny>::type create_p2l() const
	{
		return create_p2x_integral(m_fmm.template create_p2l<Ny>(),
			m_far_field_quadrature_order,
			typename type2tag<trial_field_t>::type());
	}

	/// \brief factory function for the M2P operator
	/// \tparam Nx the order of receiver differentiation
	template <int Nx>
	typename m2p_type<Nx>::type create_m2p() const
	{
		return create_x2p_integral(m_fmm.template create_m2p<Nx>(),
			m_far_field_quadrature_order,
			typename type2tag<test_field_t>::type());
	}

	/// \brief factory function for the L2P operator
	/// \tparam Nx the order of receiver differentiation
	template <int Nx>
	typename l2p_type<Nx>::type create_l2p() const
	{
		return create_x2p_integral(m_fmm.template create_l2p<Nx>(),
			m_far_field_quadrature_order,
			typename type2tag<test_field_t>::type());
	}
	
	/// \brief factory function for the P2P operator
	/// \tparam Nx the order of receiver differentiation
	/// \tparam Ny the order of source differentiation
	template <int Nx, int Ny>
	typename p2p_type<Nx, Ny>::type create_p2p() const
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
