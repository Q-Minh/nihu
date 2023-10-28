/** \file spectral_interpolate.h
 * \brief declaration of class spectral_interpolate
 */
#ifndef FMM_SPECTRAL_INTERPOLATE_H_INCLUDED
#define FMM_SPECTRAL_INTERPOLATE_H_INCLUDED

#include <fftw3.h>
#include <Eigen/Core>
#include <complex>

namespace NiHu
{
namespace fmm
{

/** \brief class performing spectral interpolation */
class spectral_interpolate
{
public:
	/** \brief the spectrum type */
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;

	/** \brief destructor */
	~spectral_interpolate();
	
	/** \brief constructor */
	spectral_interpolate(size_t Nto = 0, size_t Nfrom = 0);
	
	/** \brief copy constructor */
	spectral_interpolate(spectral_interpolate const &rhs);
	
	/** \brief move constructor */
	spectral_interpolate(spectral_interpolate &&rhs);

	/** \brief assign operator */
	spectral_interpolate &operator=(spectral_interpolate const &rhs);
	
	/** \brief move assign operator */
	spectral_interpolate &operator=(spectral_interpolate &&rhs);

	/** \brief interpolate a spectrum
	 * \param [in] from the spectrum to interpolate
	 * \return the interpolated spectrum
	 */
	cvector_t const &interpolate(cvector_t const &from) const;

private:
	void init();
	void destroy();

private:
	size_t m_Nto;
	size_t m_Nfrom;
	fftw_plan m_dft_plan;
	fftw_plan m_idft_plan;
	mutable cvector_t m_in_spec;
	mutable cvector_t m_in;
	mutable cvector_t m_out;
	mutable cvector_t m_out_spec;
};

} // end of namesapce
} // namespace NiHu

#endif // FMM_SPECTRAL_INTERPOLATE_H_INCLUDED
