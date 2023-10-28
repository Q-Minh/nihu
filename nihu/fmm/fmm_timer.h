/** 
 * @file fmm_timer.h
 * @brief Interface of the class @ref NiHu::fmm::fmm_timer
 * @ingroup fmm_util
 */
#ifndef FMM_TIMER_HPP_INCLUDED
#define FMM_TIMER_HPP_INCLUDED

#include "fmm_operator.hpp"
#include "nihu/util/timer.h"

#include <array>
#include <iostream>
#include <vector>

namespace NiHu
{
namespace fmm
{

/** \brief class to store fmm timing data */
class fmm_timer
{
	typedef NiHu::cpu_time timer_t;
	
public:
	enum {
		/** \brief index of M2M operation */
		M2M = op_tags::tag2idx(op_tags::m2m()),
		/** \brief index of L2L operation */
		L2L = op_tags::tag2idx(op_tags::l2l()),
		/** \brief index of M2L operation */
		M2L = op_tags::tag2idx(op_tags::m2l()),
		/** \brief index of P2M operation */
		P2M = op_tags::tag2idx(op_tags::p2m()),
		/** \brief index of P2L operation */
		P2L = op_tags::tag2idx(op_tags::p2l()),
		/** \brief index of L2P operation */
		L2P = op_tags::tag2idx(op_tags::l2p()),
		/** \brief index of M2P operation */
		M2P = op_tags::tag2idx(op_tags::m2p()),
		/** \brief index of P2P operation */
		P2P = op_tags::tag2idx(op_tags::p2p()),
		/** \brief number of time indices */
		NUM_TIME_INDICES = op_tags::num_tags(),
	} time_index_t;

private:
	std::vector<std::array<long long, NUM_TIME_INDICES> > m_times;
	timer_t::time_point_t m_t0;

public:
	/** \brief constructor
	 * \param [in] nLevel number of levels in the tree
	 */
	fmm_timer(size_t nLevel)
		: m_times(nLevel)
	{
	}

	/** \brief reset timer */
	void reset(void)
	{
		for (auto &it : m_times)
			it[M2M] = it[L2L] = it[M2L] = it[P2M] = it[P2L] = it[L2P] = it[M2P] = it[P2P] = 0;
	}

	/** \brief start timer */
	typename timer_t::time_point_t tic(void)
	{
		return m_t0 = timer_t::tic();
	}

	/** 
	 * \brief stop timer at a given level and operation type
	 * \param [in] level the level in the tree
	 * \param [in] type the operation type that was timed
	 */
	void toc(size_t level, int type)
	{
		m_times[level][type] += (long long) (timer_t::toc(m_t0) * 1e6);
	}

	/** 
	 * \brief inserter into output stream
	 * \param [in,out] os the output stream
	 */
	std::ostream &print(std::ostream &os = std::cout) const;

	std::vector<std::array<long long, NUM_TIME_INDICES> > const & get_times(void) const
	{
		return m_times;
	}
};

} // namespace fmm
} // namesapce NiHu

#endif // FMM_TIMER_HPP_INCLUDED
