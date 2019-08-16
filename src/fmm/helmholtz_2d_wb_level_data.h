/** \file helmholtz_2d_wb_level_data.h
 * \brief declaration of class helmholtz_2d_wb_level_data
 */
#ifndef FMM_HELMHOLTZ_2D_WB_LEVEL_DATA_HPP_INCLUDED
#define FMM_HELMHOLTZ_2D_WB_LEVEL_DATA_HPP_INCLUDED

#include <fftw3.h>
#include <Eigen/Dense>
#include <complex>
#include <vector>

#include "spectral_interpolate.h"

namespace NiHu
{
namespace fmm
{

/** \brief class containing the level data of the helmholtz 2d wide band fmm */
class helmholtz_2d_wb_level_data
{
public:
	/** \brief complex column vector type */
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;

	/** \brief constructor */
	helmholtz_2d_wb_level_data();
	/** \brief destructor */
	~helmholtz_2d_wb_level_data();
	/** \brief copy constructor */
	helmholtz_2d_wb_level_data(helmholtz_2d_wb_level_data const &);
	/** \brief assignment operator */
	helmholtz_2d_wb_level_data &operator=(helmholtz_2d_wb_level_data const &) = delete;

	/** \brief initialize the data for a given frequency
	 * \param [in] drel the relative cluster size (in wave lengths)
	 */
	void init(double drel);

	/** \brief set the expansion length of the level
	 * \param [in] expansion_length the expansion length
	 */
	void set_expansion_length(size_t expansion_length);
	
	/** \brief return the expansion length */
	size_t get_expansion_length() const;
	
	/** \brief return the size of the multipole / local data */
	size_t get_data_size() const;

	/** \brief set if the level is in high frequency */
	void set_high_freq(bool hf);
	
	/** \brief return true if the level is in hgh frequency */
	bool is_high_freq() const;

	/** \brief compute upward interpolating matrix
	 \param [in] Nfrom the expansion length in the lower cluster
	 */
	void set_interp_up(size_t  Nfrom);
	
	/** \brief compute the downward interpolating matrix
	 \param [in] Nfrom the expansion length in the upper cluster
	 */
	void set_interp_dn(size_t Nfrom);

	/** \brief interpolate multipole contribution up to this level
	 * \param [in] multi the multipole coefficient in the lower level
	 * \return the multipole coefficient on this level
	 */
	cvector_t const &interp_up(cvector_t const &multi) const;

	/** \brief interpolate local data down to this level
	 * \param [in] local the local coefficient in the higher level
	 * \return the local coefficient on this level
	 */
	cvector_t const &interp_dn(cvector_t const &local) const;

	/** \brief perform dft on this level
	 * \param [in] multi the multipole coefficient vector in coeff domain
	 * \return the multipole coefficient vector in spectral domain
	 */
	cvector_t const &dft(cvector_t const &multi) const;
	
	/** \brief perform inverse dft on this level
	 * \param [in] local the local coefficient vector in spectral domain
	 * \return the local coefficient vector in coefficient domain
	 */
	cvector_t const &idft(cvector_t const &local) const;

	void set_num_threads(size_t n);

private:
	size_t m_expansion_length;
	bool m_high_freq;

	std::vector<spectral_interpolate> m_interp_ups;	// should be thread private
	std::vector<spectral_interpolate> m_interp_dns;	// should be thread private

	fftw_plan m_dft_plan;
	mutable std::vector<cvector_t> m_multi_paddeds;	// should be thread private
	mutable std::vector<cvector_t> m_multi_spectrums;	// should be thread private

	fftw_plan m_idft_plan;
	mutable std::vector <cvector_t> m_local_paddeds;	// should be thread private
	mutable std::vector <cvector_t> m_locals;			// should be thread private
};

} // end of namespace fmm
} // namespace NiHu

#endif
