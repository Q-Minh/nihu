/** \file helmholtz_2d_wb_x2x_matrix.cpp
 * \brief implementation of the m2m l2l and m2l operators for the helmholtz 2d wide band fmm
 */
#include "helmholtz_2d_wb_x2x_matrix.h"

#include <iostream>

namespace NiHu
{
namespace fmm
{

helmholtz_2d_wb_m2m_matrix::helmholtz_2d_wb_m2m_matrix()
	: m_level_data_from(nullptr)
	, m_level_data_to(nullptr)
{
}


helmholtz_2d_wb_m2m_matrix::helmholtz_2d_wb_m2m_matrix(
	helmholtz_2d_wb_level_data const &level_to,
	helmholtz_2d_wb_level_data const &level_from,
	cvector_t diag_coeffs)
	: m_level_data_to(&level_to)
	, m_level_data_from(&level_from)
	, m_cm(level_to.get_expansion_length(), level_from.get_expansion_length(), diag_coeffs)
{
	if (m_level_data_from->get_high_freq() && m_level_data_to->get_high_freq())
		m_transfer = m_level_data_to->dft(diag_coeffs);
}

helmholtz_2d_wb_m2m_matrix::cvector_t
helmholtz_2d_wb_m2m_matrix::operator*(cvector_t const &rhs) const
{
	// low to low
	if (!m_level_data_from->get_high_freq() && !m_level_data_to->get_high_freq())
		return m_cm * rhs;

	// low to high
	if (!m_level_data_from->get_high_freq() && m_level_data_to->get_high_freq())
	{
		cvector_t multi = m_cm * rhs;
		cvector_t spect = m_level_data_to->dft(multi);
		return spect;
	}

	// high to high
	return m_level_data_to->interp_up(rhs).array() * m_transfer.array();
}


helmholtz_2d_wb_l2l_matrix::helmholtz_2d_wb_l2l_matrix()
	: m_level_data_from(nullptr)
	, m_level_data_to(nullptr)
{
}

helmholtz_2d_wb_l2l_matrix::helmholtz_2d_wb_l2l_matrix(helmholtz_2d_wb_level_data const &level_to,
	helmholtz_2d_wb_level_data const &level_from,
	cvector_t diag_coeffs)
	: m_level_data_to(&level_to)
	, m_level_data_from(&level_from)
	, m_cm(level_to.get_expansion_length(), level_from.get_expansion_length(), diag_coeffs)
{
	if (m_level_data_from->get_high_freq() && m_level_data_to->get_high_freq())
		m_transfer = m_level_data_from->dft(diag_coeffs);
}

helmholtz_2d_wb_l2l_matrix::cvector_t
helmholtz_2d_wb_l2l_matrix::operator*(cvector_t const &rhs) const
{
	if (!m_level_data_from->get_high_freq() && !m_level_data_to->get_high_freq())
		return m_cm * rhs;

	// high to low
	if (m_level_data_from->get_high_freq() && !m_level_data_to->get_high_freq())
	{
		cvector_t local_padded = m_level_data_from->idft(rhs);
		cvector_t spect = m_level_data_from->dft(local_padded);
		cvector_t local = m_cm * local_padded;
		return local;
	}

	// L2L high to high
	return m_level_data_to->interp_dn(m_transfer.array() * rhs.array());
}


helmholtz_2d_wb_m2l_matrix::helmholtz_2d_wb_m2l_matrix()
	: m_level_data(nullptr)
{
}


helmholtz_2d_wb_m2l_matrix::helmholtz_2d_wb_m2l_matrix(helmholtz_2d_wb_level_data const &level,
	cvector_t diag_coeffs)
	: m_level_data(&level)
	, m_cm(m_level_data->get_expansion_length(), m_level_data->get_expansion_length(), diag_coeffs)
{
	if (m_level_data->get_high_freq())
		m_transfer = m_level_data->dft(diag_coeffs);
}

helmholtz_2d_wb_m2l_matrix::cvector_t
helmholtz_2d_wb_m2l_matrix::operator*(cvector_t const &rhs) const
{
	if (!m_level_data->get_high_freq())
		return m_cm * rhs;

	return rhs.array() * m_transfer.array();
}

} // end of namespace fmm
} // namespace NiHu