/** 
 * \file operator_with_wave_number.hpp
 * \brief implementation of class \ref NiHu::fmm::operator_with_wave_number
 * \ingroup fmm_ops
 */

#ifndef OPERATOR_WITH_WAVE_NUMBER_HPP_INCLUDED
#define OPERATOR_WITH_WAVE_NUMBER_HPP_INCLUDED

namespace NiHu
{
namespace fmm
{
/// \brief class storing a wave number
/// \tparam WaveNumber the wave number type
template <class WaveNumber>
class operator_with_wave_number
{
public:
	typedef WaveNumber wave_number_t;

	operator_with_wave_number(wave_number_t const &wave_number)
		: m_wave_number(wave_number)
	{
	}

	wave_number_t const &get_wave_number() const
	{
		return this->m_wave_number;
	}

private:
	wave_number_t const m_wave_number;
};
	
} // end of namespace fmm
} // end namespace NiHu

#endif /* OPERATOR_WITH_WAVE_NUMBER_HPP_INCLUDED */
