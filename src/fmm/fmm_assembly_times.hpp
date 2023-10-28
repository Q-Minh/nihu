#ifndef FMM_ASSEMBLY_TIMES_HPP_INCLUDED
#define FMM_ASSEMBLY_TIMES_HPP_INCLUDED

#include "fmm_timer.h"

/**
 * @brief Helper class for storing FMM assembly times 
 */
class fmm_assembly_times
{
	using timer = NiHu::fmm::fmm_timer;
public:
	/** 
	 * @brief Fill the assembly time from an operator collection
	 * @tparam Collection The operator collection class
	 * @param[in] coll Collection instance
	 */
	template <class Collection>
	void fill_times(Collection const & coll)
	{
		m_times[timer::M2M] = coll.get(NiHu::fmm::op_tags::idx2tag<timer::M2M>::type()).get_assembly_time();
		m_times[timer::L2L] = coll.get(NiHu::fmm::op_tags::idx2tag<timer::L2L>::type()).get_assembly_time();
		m_times[timer::M2L] = coll.get(NiHu::fmm::op_tags::idx2tag<timer::M2L>::type()).get_assembly_time();
		
		m_times[timer::P2M] = coll.get(NiHu::fmm::op_tags::idx2tag<timer::P2M>::type()).get_assembly_time();
		m_times[timer::P2L] = coll.get(NiHu::fmm::op_tags::idx2tag<timer::P2L>::type()).get_assembly_time();
		m_times[timer::L2P] = coll.get(NiHu::fmm::op_tags::idx2tag<timer::L2P>::type()).get_assembly_time();
		m_times[timer::M2P] = coll.get(NiHu::fmm::op_tags::idx2tag<timer::M2P>::type()).get_assembly_time();
		
		m_times[timer::P2P] = coll.get(NiHu::fmm::op_tags::idx2tag<timer::P2P>::type()).get_assembly_time();
	}
	
	/**
	 * @brief Returns the assembly time for the given operation index
	 * @return Assembly time in microsecond units
	 */
	size_t const get_time(size_t idx) const
	{
		return m_times[idx];
	}
	
private:
	size_t m_times[timer::NUM_TIME_INDICES];
};

#endif /* FMM_ASSEMBLY_TIMES_HPP_INCLUDED */
