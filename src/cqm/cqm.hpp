/** \file cqm.hpp
 * \brief A generic implementation of the Convolution Quadrature Method (CQM)
 * \author Peter Fiala fiala@hit.bme.hu
 */

#ifndef CQM_HPP_INCLUDED
#define CQM_HPP_INCLUDED

#include "../util/math_constants.hpp"

#include <complex>
#include <vector>
#include <algorithm>	// for std::transform
#include <fftw3.h>

/** \brief Backward-differentiation multistep method ratio
 * \tparam Order the order of the backward differentiation
 */
template <unsigned Order>
struct BDF;

/** \brief Specialisation of BDF for 1st order */
template <>
struct BDF<1>
{
	/** \brief return ratio in z domain
	 * \tparam T the scalar type
	 * \return the BDF ratio in z domain
	 */
	template <class T>
	static T gamma(T const &z) { return 1. - z; }
};

/** \brief Specialisation of BDF for 2nd order */
template <>
struct BDF<2>
{
	/** \brief return ratio in z domain
	 * \tparam T the scalar type
	 * \return the BDF ratio in z domain
	 */
	template <class T>
	static T gamma(T const &z) { return 1.5 + z*(-2. + .5*z); }
};



/** \brief Convolution Quadrature implementation
 * \tparam KernelResult the result type of the kernel in the Laplace domain
 * \tparam Order the BDF order
 */
template <unsigned Order, class ExcTime, class ExcFreq, class RespTime, class RespFreq>
class CQM
{
public:
	/** \brief template argument as nested constant */
	enum { BDF_order = Order };

	/** \brief constructor
	 * \param N the convolution length
	 * \param dt the time step
	 * \param delta accuracy of kernel evaluation in Laplace domain
	 */
	CQM(int N, double dt, double delta) :
		m_N(N),
		m_dt(dt),
		m_delta(delta),
		m_rho(std::exp(std::log(m_delta) / (2.*m_N))),
		m_svec(m_N / 2 + 1),
		m_scale(m_N, 0.),
		m_scaled_excitation(m_N),
		m_laplace_excitation(m_N / 2 + 1),
		m_laplace_response(m_N / 2 + 1),
		m_time_response(m_N)
	{
		double dphi = 2.*M_PI / m_N;
		double phi = 0.;
		// conjugated because fftw r2c computes the forward transform with negative sign
		for (unsigned k = 0; k < m_N / 2 + 1; ++k, phi += dphi)
			m_svec[k] = BDF<BDF_order>::gamma(m_rho * std::complex<double>(cos(phi), -sin(phi)) ) / m_dt;

		int n[] = {m_N};
		int rank = 1;
		int xstride = (double *)(m_scaled_excitation.data()+1) - (double *)m_scaled_excitation.data();
		int rstride = (double *)(m_time_response.data()+1) - (double *)m_time_response.data();
		int const xhowmany = xstride;	// this is the default only
		int const rhowmany = rstride;
		int xdist = 1, rdist = 1;

		// design fft plan for ifft
		m_ifftw_plan = fftw_plan_many_dft_r2c(
			rank, n, xhowmany,
			reinterpret_cast<double *>(m_scaled_excitation.data()), NULL,
			xstride, xdist,
			reinterpret_cast<fftw_complex *>(m_laplace_excitation.data()), NULL,
			xstride, xdist,
			FFTW_MEASURE);

		m_fftw_plan = fftw_plan_many_dft_c2r(
			rank, n, rhowmany,
			reinterpret_cast<fftw_complex *>(m_laplace_response.data()), NULL,
			rstride, rdist,
			reinterpret_cast<double *>(m_time_response.data()), NULL,
			rstride, rdist,
			FFTW_MEASURE);

		// compute samples of scaling function
		double r = 1.;
		for (auto it = m_scale.begin(); it != m_scale.end(); ++it, r *= m_rho)
			*it = r;
	}

	/** \brief destructor */
	~CQM()
	{
		fftw_destroy_plan(m_fftw_plan);
		fftw_destroy_plan(m_ifftw_plan);
	}

	/** \brief return convolution length */
	unsigned get_N(void) const { return m_N; }
	
	/** \brief return time step */
	double get_dt(void) const { return m_dt; }
	
	/** \brief return kernel accuracy */
	double get_delta(void) const { return m_delta; }
	
	/** \brief return radius of integration on the complex plane */
	double get_rho(void) const { return m_rho; }
	
	/** \brief return vector of complex frequencies */
	std::vector<std::complex<double> > const &get_svec(void) const { return m_svec; }
	
	/** \brief return discrete impulse response result */
	std::vector<std::complex<RespFreq> > const &get_laplace_response(void) const { return m_laplace_response; }
	
	/** \brief return discrete impulse response result */
	std::vector<RespTime> const &get_time_response(void) const { return m_time_response; }

	template <class ExcIt, class Kernel>
	void eval(ExcIt begin, ExcIt end, Kernel const &kernel)
	{
		// scale-copy excitation to scaled excitation
		for (int i = 0; i < m_N; ++i)
			m_scaled_excitation[i] = m_scale[i] * *(begin+i);
		// ifft of scaled excitation
		fftw_execute(m_ifftw_plan);
		// evaluate kernel and multiply with excitation
		for (int i = 0; i < m_N / 2 + 1; ++i)
			m_laplace_response[i] = kernel(m_svec[i]).eval() * m_laplace_excitation[i];
		// fft of result
		fftw_execute(m_fftw_plan);
		// rescale in time domain
		for (int i = 0; i < m_N; ++i)
			m_time_response[i] /= (m_N * m_scale[i]);
	}

private:
	int m_N;				/** \brief convolution depth */
	double m_dt;			/** \brief time step */
	double m_delta;			/** \brief kernel accuracy */
	double m_rho;			/** \brief radius of integration on s plane */
	fftw_plan m_fftw_plan;	/** \brief DFT plan */
	fftw_plan m_ifftw_plan;	/** \brief DFT plan */

	std::vector<std::complex<double> > m_svec;		/** \brief complex frequency samples */
	std::vector<double> m_scale;					/** \brief scale function */
	std::vector<ExcTime> m_scaled_excitation;		/** \brief scaled excitation samples */
	std::vector<ExcFreq> m_laplace_excitation;		/** \brief s-domain excitation */
	std::vector<RespFreq> m_laplace_response;		/** \brief s domain discrete response */
	std::vector<RespTime> m_time_response;			/** \brief time domain discrete response */
};

#endif // CQM_HPP_INCLUDED

