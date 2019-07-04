#include "fmm_timer.hpp"

#include <iomanip>
#include <sstream>

namespace NiHu
{
namespace fmm
{

std::ostream &fmm_timer::print(std::ostream &os) const
{
	// save state of os
	std::stringstream state;
	state.copyfmt(os);

	os << std::fixed << std::setprecision(1);

	size_t nLevel = m_times.size();
	for (size_t iLevel = 0; iLevel < nLevel; ++iLevel)
	{
		os << "Level #" << std::setw(2) << iLevel << ": "
			<< "M2M (" << std::setw(10) << m_times[iLevel][M2M] << "), "
			<< "M2L (" << std::setw(10) << m_times[iLevel][M2L] << "), "
			<< "L2L (" << std::setw(10) << m_times[iLevel][L2L] << ")" << std::endl;
	}

	os << "P2P (" << m_times[0][P2P] << "), "
		<< "P2M (" << m_times[0][P2M] << "), "
		<< "P2L (" << m_times[0][P2L] << "), "
		<< "L2P (" << m_times[0][L2P] << "), "
		<< "M2P (" << m_times[0][M2P] << ")" << std::endl;

	// restore state of os
	os.copyfmt(state);

	return os;
}

} // namespace fmm
} // namespace NiHu
