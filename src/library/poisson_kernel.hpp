/**
 * \file poisson_kernel.hpp
 * \ingroup library
 * \brief implementation of kernels of the Poisson equation \f$\nabla^2 p = 0\f$
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef POISSON_KERNEL_HPP_INCLUDED
#define POISSON_KERNEL_HPP_INCLUDED

#include "../bem/kernel.hpp"
#include "../bem/gaussian_quadrature.hpp"

/** \brief 3D Poisson kernel \f$1/4\pi r\f$ */
class poisson_G_kernel;

/** \brief traits of the poisson G kernel */
template<>
struct kernel_traits<poisson_G_kernel>
{
	/** \brief kernel test input type */
	typedef location<space_3d> test_input_t;
	/** \brief kernel trial input type */
	typedef location<space_3d> trial_input_t;
	/** \brief kernel result type */
	typedef space_3d::scalar_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 1;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
};

/** \brief 3D Poisson kernel \f$1/4\pi r\f$ */
class poisson_G_kernel :
	public kernel_base<poisson_G_kernel>
{
public:
	/** \brief evaluate kernel between test and trial positions
	* \param [in] x the test input
	* \param [in] y the trial input
	* \return the kernel value K(x,y)
	*/
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		scalar_t r = (y.get_x() - x.get_x()).norm();
		return 1.0 / r / (4.0 * M_PI);
	}
	
	/** \brief estimate kernel's polynomial complexity
	* \param [in] x the test input
	* \param [in] y the trial input
	* \param [in] s estimate linear size of one of the elements
	* \return the kernel value K(x,y)
	*/
	unsigned estimate_complexity(
		test_input_t const &x,
		trial_input_t const &y,
		scalar_t s) const
	{
		/** \todo hard coding of this 3 is sick */
		return 3;
	}
};


/** \brief 3D Poisson kernel \f$-1/4\pi r^2 \cdot dr/dn \f$ */
class poisson_H_kernel;

/** \brief traits of the Poisson H kernel */
template<>
struct kernel_traits<poisson_H_kernel>
{
	/** \brief kernel test input type */
	typedef location<space_3d> test_input_t;
	/** \brief kernel trial input type */
	typedef location_with_normal<space_3d> trial_input_t;
	/** \brief kernel result type */
	typedef space_3d::scalar_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 2;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
};

/** \brief 3D Poisson kernel \f$ -1/4\pi r^2 \cdot dr/dn \f$ */
class poisson_H_kernel :
	public kernel_base<poisson_H_kernel>
{
public:
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

		result_t m_result = 1.0 / r / (4.0 * M_PI);
		scalar_t rdn = rvec.dot(y.get_unit_normal());
		m_result *= (-rdn / r2);

		return m_result;
	}

	/** \brief estimate kernel's polynomial complexity
	* \param [in] x the test input
	* \param [in] y the trial input
	* \param [in] s estimate linear size of one of the elements
	* \return the kernel value K(x,y)
	*/
	unsigned estimate_complexity(
		test_input_t const &x,
		trial_input_t const &y,
		scalar_t s) const
	{
		/** \todo hard coding of this 5 is sick */
		return 5;
	}
};


/** \brief 3D Poisson kernel \f$ 1/4\pi r \left\{1, -1/r \cdot dr/dn\right\} \f$ */
class poisson_GH_kernel;

/** \brief traits of the double Poisson kernel */
template<>
struct kernel_traits<poisson_GH_kernel>
{
	/** \brief kernel test input type */
	typedef location<space_3d> test_input_t;
	/** \brief kernel trial input type */
	typedef location_with_normal<space_3d> trial_input_t;
	/** \brief kernel result type */
	typedef couple<space_3d::scalar_t> result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 1;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
};

/** \brief 3D Poisson kernel \f$ 1/4\pi r \left\{1, -1/r \cdot dr/dn\right\} \f$ */
class poisson_GH_kernel :
	public kernel_base<poisson_GH_kernel>
{
public:
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

		result_t m_result(1.0 / r / (4.0 * M_PI));
		scalar_t rdn = rvec.dot(y.get_unit_normal());
		m_result.second() = m_result.first() * (-rdn / r2);

		return m_result;
	}
	
	/** \brief estimate kernel's polynomial complexity
	* \param [in] x the test input
	* \param [in] y the trial input
	* \param [in] s estimate linear size of one of the elements
	* \return the kernel value K(x,y)
	*/
	unsigned estimate_complexity(
		test_input_t const &x,
		trial_input_t const &y,
		scalar_t s) const
	{
		/** \todo hard coding of this 5 is sick */
		return 5;
	}
};

#endif // POISSON_KERNEL_HPP_INCLUDED


