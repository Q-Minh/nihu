/** \file fmm_timer.hpp
 * \brief declaration of class fmm_timer
 */
#ifndef FMM_TIMER_HPP_INCLUDED
#define FMM_TIMER_HPP_INCLUDED

#include "util/timer.h"

#include <vector>
#include <array>
#include <iostream>

namespace NiHu
{
namespace fmm
{

/** \brief class to store fmm timing data */
class fmm_timer
{
	typedef NiHu::cpu_time timer_t;
	
public:
	/** \brief index of M2M operation */
	static int const M2M = 0;
	/** \brief index of L2L operation */
	static int const L2L = 1;
	/** \brief index of M2L operation */
	static int const M2L = 2;
	/** \brief index of P2M operation */
	static int const P2M = 3;
	/** \brief index of P2L operation */
	static int const P2L = 4;
	/** \brief index of L2P operation */
	static int const L2P = 5;
	/** \brief index of M2P operation */
	static int const M2P = 6;
	/** \brief index of P2P operation */
	static int const P2P = 7;

private:
	std::vector<std::array<long long, 8> > m_times;
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

	/** \brief stop timer at a given level and operation type
	 * \param [in] level the level in the tree
	 * \param [in] type the operation type that was timed
	 */
	void toc(size_t level, int type)
	{
		m_times[level][type] += (long long) (timer_t::toc(m_t0) * 1e6);
	}

	/** \brief inserter into output stream
	 * \param [in, out] os the output stream
	 */
	std::ostream &print(std::ostream &os = std::cout) const;

	std::vector<std::array<long long, 8> > const & get_times(void) const
	{
		return m_times;
	}
};

} // namespace fmm
} // namesapce NiHu

#endif // FMM_TIMER_HPP_INCLUDED
