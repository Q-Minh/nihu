/** 
 * \file helmholtz_2d_wb_x2x_matrix.h
 * \brief declaration of m2m l2l and m2l operations for the 2d wide band helmholtz fmm
 * \ingroup hmm_helmholtz_2d_wb
 */

#ifndef FMM_HELMHOLTZ_2D_WB_X2X_MATRIX_H_INCLUDED
#define FMM_HELMHOLTZ_2D_WB_X2X_MATRIX_H_INCLUDED

#include "convolution_matrix.hpp"
#include "helmholtz_2d_wb_level_data.h"

namespace NiHu
{
namespace fmm
{

/** \brief m2m matrix of the wide band 2d helmholtz fmm */
class helmholtz_2d_wb_m2m_matrix
{
public:
	/** \brief complex column vector (the multipole type) */
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;

	/** \brief costructor */
	helmholtz_2d_wb_m2m_matrix();

	/** \brief constructor
	 * \param[in] level_data_to pointer to the receiver level data
	 * \param[in] level_data_from pointer to the source level data
	 * \param[in] diag_coeffs diagonal coefficients of the convolution
	 */
	helmholtz_2d_wb_m2m_matrix(
		helmholtz_2d_wb_level_data const &level_data_to,
		helmholtz_2d_wb_level_data const &level_data_from,
		cvector_t diag_coeffs);

	/** \brief multiply the matrix with a multipole coefficient from the right
	 * \param [in] rhs the right hand side multipole coefficient
	 * \return the resulting multipole coefficient
	 */
	cvector_t operator*(cvector_t const &rhs) const;

private:
	helmholtz_2d_wb_level_data const *m_level_data_to;
	helmholtz_2d_wb_level_data const *m_level_data_from;
	convolution_matrix<std::complex<double> > m_cm;
	cvector_t m_transfer;
};


/** \brief l2l matrix of the wide band 2d helmholtz fmm */
class helmholtz_2d_wb_l2l_matrix
{
public:
	/** \brief complex column vector (the local type) */
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;

	/** \brief costructor */
	helmholtz_2d_wb_l2l_matrix();

	/** \brief constructor
	 * \param[in] level_data_to pointer to the receiver level data
	 * \param[in] level_data_from pointer to the source level data
	 * \param[in] diag_coeffs diagonal coefficients of the convolution
	 */
	helmholtz_2d_wb_l2l_matrix(
		helmholtz_2d_wb_level_data const &level_data_to,
		helmholtz_2d_wb_level_data const &level_data_from,
		cvector_t diag_coeffs);

	/** \brief multiply the matrix with a local coefficient from the right
	 * \param [in] rhs the right hand side local coefficient
	 * \return the resulting local coefficient
	 */
	cvector_t operator*(cvector_t const &rhs) const;

private:
	helmholtz_2d_wb_level_data const *m_level_data_to;
	helmholtz_2d_wb_level_data const *m_level_data_from;
	convolution_matrix<std::complex<double> > m_cm;
	cvector_t m_transfer;
};


/** \brief m2l matrix of the wide band 2d helmholtz fmm */
class helmholtz_2d_wb_m2l_matrix
{
public:
	/** \brief complex column vector (the multipole and local type) */
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;

	/** \brief costructor */
	helmholtz_2d_wb_m2l_matrix();

	/** \brief constructor
	 * \param[in] level_data pointer to the level data
	 * \param[in] diag_coeffs diagonal coefficients of the convolution
	 */
	helmholtz_2d_wb_m2l_matrix(
		helmholtz_2d_wb_level_data const &level_data,
		cvector_t diag_coeffs);

	/** \brief multiply the matrix with a multipole coefficient from the right
	 * \param [in] rhs the right hand side multipole coefficient
	 * \return the resulting local coefficient
	 */
	cvector_t operator*(cvector_t const &rhs) const;

private:
	helmholtz_2d_wb_level_data const *m_level_data;
	convolution_matrix<std::complex<double> > m_cm;
	cvector_t m_transfer;
};

} // end of namespace fmm
} // namespace NiHu
#endif
