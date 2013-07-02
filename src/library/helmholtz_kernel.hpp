/**
 * \file helmholtz_kernel.hpp
 * \ingroup library
 * \brief implementation of kernels of the Helmholtz equation \f$\nabla^2 p + k^2 p = 0\f$
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef HELMHOLTZ_KERNEL_HPP_INCLUDED
#define HELMHOLTZ_KERNEL_HPP_INCLUDED

#include "../bem/kernel.hpp"
#include "../bem/gaussian_quadrature.hpp"

#include "location_and_normal.hpp"
#include "reciprocal_distance_kernel.hpp"

/** \brief Common helper base for a kernel with a wave number member
* \tparam wave_number_type the type of the wave number
*/
template <class wave_number_type>
class kernel_with_wave_number
{
public:
	/** \brief set wave number to a value
	* \param [in] k the new wave number
	*/
	void set_wave_number(wave_number_type const &k)
	{
		m_k = k;
	}
	
	/** \brief get the wave number
	* \return the new wave number
	*/
	wave_number_type const &get_wave_number(void) const
	{
		return m_k;
	}

private:
	/** \brief the stored wave number */
	wave_number_type m_k;
};


/** \brief 3D Helmholtz kernel \f$\exp(-ikr)/4\pi r\f$ */
class helmholtz_G_kernel;

/** \brief traits of the helmholtz G kernel */
template<>
struct kernel_traits<helmholtz_G_kernel>
{
	/** \brief kernel test input type */
	typedef location<space_3d> test_input_t;
	/** \brief kernel trial input type */
	typedef location<space_3d> trial_input_t;
	/** \brief kernel result type */
	typedef std::complex<space_3d::scalar_t> result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 1;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
};

/** \brief 3D Helmholtz kernel \f$\exp(-ikr)/4\pi r\f$ */
class helmholtz_G_kernel :
	public kernel_base<helmholtz_G_kernel>,
	public kernel_with_wave_number<std::complex<kernel_base<helmholtz_G_kernel>::scalar_t> >,
	public reciprocal_distance_kernel<helmholtz_G_kernel>
{
public:
	typedef kernel_base<helmholtz_G_kernel> base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::scalar_t scalar_t;

	/** \brief evaluate kernel between test and trial positions
	* \param [in] x the test input
	* \param [in] y the trial input
	* \return the kernel value K(x,y)
	*/
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		scalar_t r = (y.get_x() - x.get_x()).norm();
		return std::exp(-std::complex<scalar_t>(0.0,1.0)*get_wave_number()*r) / r / (4.0 * M_PI);
	}
};


/** \brief 3D Helmholtz kernel \f$ \exp(-ikr)/4\pi r \left(-(1+ikr)/r\right) \cdot dr/dn \f$ */
class helmholtz_H_kernel;

/** \brief traits of the Helmholtz H kernel */
template<>
struct kernel_traits<helmholtz_H_kernel>
{
	/** \brief kernel test input type */
	typedef location<space_3d> test_input_t;
	/** \brief kernel trial input type */
	typedef location_with_normal<space_3d> trial_input_t;
	/** \brief kernel result type */
	typedef std::complex<space_3d::scalar_t> result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 2;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
};

/** \brief 3D Helmholtz kernel \f$ \exp(-ikr)/4\pi r \left(-(1+ikr)/r\right) \cdot dr/dn \f$ */
class helmholtz_H_kernel :
	public kernel_base<helmholtz_H_kernel>,
	public kernel_with_wave_number<std::complex<kernel_base<helmholtz_G_kernel>::scalar_t> >,
	public reciprocal_distance_kernel<helmholtz_H_kernel>
{
public:
	typedef kernel_base<helmholtz_H_kernel> base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::scalar_t scalar_t;

	/** \brief evaluate kernel between test and trial positions
	* \param [in] x the test input
	* \param [in] y the trial input
	* \return the kernel value K(x,y)
	*/
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		x_t rvec = y.get_x() - x.get_x();
		scalar_t r2 = rvec.squaredNorm();
		scalar_t r = sqrt(r2);
		std::complex<scalar_t> ikr(std::complex<scalar_t>(0.,1.)*get_wave_number()*r);

		result_t m_result = std::exp(-ikr) / r / (4.0 * M_PI);
		scalar_t rdn = rvec.dot(y.get_unit_normal());
		m_result *= (1.0 + ikr) * (-rdn / r2);

		return m_result;
	}
};


/** \brief 3D Helmholtz kernel \f$ \exp(-ikr)/4\pi r \left\{1, -(1+ikr)/r \cdot dr/dn\right\} \f$ */
class helmholtz_GH_kernel;

/** \brief traits of the double helmholtz kernel */
template<>
struct kernel_traits<helmholtz_GH_kernel>
{
	/** \brief kernel test input type */
	typedef location<space_3d> test_input_t;
	/** \brief kernel trial input type */
	typedef location_with_normal<space_3d> trial_input_t;
	/** \brief kernel result type */
	typedef couple<std::complex<space_3d::scalar_t> > result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 2;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
};

/** \brief 3D Helmholtz kernel \f$ \exp(-ikr)/4\pi r \left\{1, -(1+ikr)/r \cdot dr/dn\right\} \f$ */
class helmholtz_GH_kernel :
	public kernel_base<helmholtz_GH_kernel>,
	public kernel_with_wave_number<std::complex<kernel_base<helmholtz_G_kernel>::scalar_t> >,
	public reciprocal_distance_kernel<helmholtz_GH_kernel>
{
public:
	typedef kernel_base<helmholtz_GH_kernel> base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::scalar_t scalar_t;

	/** \brief evaluate kernel between test and trial positions
	* \param [in] x the test input
	* \param [in] y the trial input
	* \return the kernel value K(x,y)
	*/
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		x_t rvec = y.get_x() - x.get_x();
		scalar_t r2 = rvec.squaredNorm();
		scalar_t r = sqrt(r2);
		std::complex<scalar_t> ikr(std::complex<scalar_t>(0.,1.)*get_wave_number()*r);

		result_t m_result(std::exp(-ikr) / r / (4.0 * M_PI));
		scalar_t rdn = rvec.dot(y.get_unit_normal());
		m_result.second() = m_result.first() * (1.0 + ikr) * (-rdn / r2);

		return m_result;
	}
};

#endif // HELMHOLTZ_KERNEL_HPP_INCLUDED


