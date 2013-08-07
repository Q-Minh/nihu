/**
 * \file helmholtz_kernel.hpp
 * \ingroup library
 * \brief implementation of kernels of the Helmholtz equation \f$ \nabla^2 p + k^2 p = 0 \f$
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef HELMHOLTZ_KERNEL_HPP_INCLUDED
#define HELMHOLTZ_KERNEL_HPP_INCLUDED

#include <cmath>
#include <complex>

#include "../bem/kernel.hpp"
#include "../bem/gaussian_quadrature.hpp"

#include "location_normal.hpp"
#include "reciprocal_distance_kernel.hpp"

#include "basic_bricks.hpp"


/** \brief a brick representing the expression \f$ i k r \f$
 * \tparam scalar scalar type of the coordinate space the distance is defined over
 */
template <class scalar>
struct ikr_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		typedef std::complex<scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel the kernel instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_t const &kernel) :
			wall(test_input, trial_input, kernel),
			m_ikr(std::complex<scalar>(0.0,1.0) * kernel.derived().get_wave_number() * wall::get_distance())
		{
		}

		/** \brief return ikr
		 * \return ikr
		 */
		std::complex<scalar> const & get_ikr(void) const
		{
			return m_ikr;
		}

	private:
		std::complex<scalar> m_ikr;
	};
};


/** \brief a brick representing a 3D Helmholtz kernel \f$ 1/4\pi r \f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct helmholtz_3d_g_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		typedef std::complex<scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel the kernel instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_t const &kernel) :
			wall(test_input, trial_input, kernel),
			m_helmholtz_g(std::exp(-1.0 * wall::get_ikr()) / wall::get_distance() / (4.0 * M_PI))
		{
		}

		/** \brief return helmholtz g kernel
		 * \return helmholtz g kernel
		 */
		std::complex<scalar> const & get_helmholtz_g(void) const
		{
			return m_helmholtz_g;
		}

		/** \brief return Helmholtz g kernel
		 * \return Helmholtz g kernel
		 */
		std::complex<scalar> const & get_result(void) const
		{
			return m_helmholtz_g;
		}

	private:
		std::complex<scalar> m_helmholtz_g;
	};
};


/** \brief combination of ::distance_vector_brick, ::distance_brick, ::ikr_brick and ::helmholtz_2d_g_brick into a wall
 * \tparam space the coordinate space the Helmholtz kernel is defined over
 */
template <class scalar>
struct helmholtz_3d_g_wall : build<
	distance_vector_brick<space<scalar, 3> >,
	distance_brick<scalar>,
	ikr_brick<scalar>,
	helmholtz_3d_g_brick<scalar>
> {};

/** \brief 3D Helmholtz kernel \f$ \exp(-ikr)/4\pi r \f$ */
class helmholtz_3d_G_kernel;

/** \brief traits of the helmholtz G kernel */
template<>
struct kernel_traits<helmholtz_3d_G_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d> >::type trial_input_t;
	/** \brief the kernel output type */
	typedef helmholtz_3d_g_wall<space_3d::scalar_t>::type output_t;
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


/** \brief 3D Helmholtz G kernel \f$ 1/4\pi r\f$ */
class helmholtz_3d_G_kernel :
	public kernel_base<helmholtz_3d_G_kernel>,
	public reciprocal_distance_kernel<helmholtz_3d_G_kernel>
{
public:
	using reciprocal_distance_kernel<helmholtz_3d_G_kernel>::estimate_complexity;

	helmholtz_3d_G_kernel(std::complex<double> wave_number) :
		m_wave_number(wave_number)
	{
	}

	std::complex<double> const &get_wave_number(void) const
	{
		return m_wave_number;
	}

private:
	std::complex<double> m_wave_number;
};

typedef helmholtz_3d_G_kernel helmholtz_3d_SLP_kernel;



/** \brief a brick representing a 3D Helmholtz H kernel \f$ \exp(-ikr)/4\pi r \left(-(1+ikr)/r\right) \cdot dr/dn \f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct helmholtz_3d_h_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		typedef std::complex<scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel the kernel instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_t const &kernel) :
			wall(test_input, trial_input, kernel),
			m_helmholtz_h(-(1.0+wall::get_ikr()) * wall::get_helmholtz_g() / wall::get_distance() * wall::get_rdny())
		{
		}

		/** \brief return helmholtz g kernel
		 * \return helmholtz g kernel
		 */
		std::complex<scalar> const & get_helmholtz_h(void) const
		{
			return m_helmholtz_h;
		}

		/** \brief return Helmholtz h kernel
		 * \return Helmholtz h kernel
		 */
		std::complex<scalar> const & get_result(void) const
		{
			return m_helmholtz_h;
		}

	private:
		std::complex<scalar> m_helmholtz_h;
	};
};


template <class scalar>
struct helmholtz_3d_h_wall : build<
	distance_vector_brick<space<scalar, 3> >,
	distance_brick<scalar>,
	ikr_brick<scalar>,
	helmholtz_3d_g_brick<scalar>,
	rdny_brick<scalar>,
	helmholtz_3d_h_brick<scalar>
> {};


/** \brief 3D Helmholtz kernel \f$ \exp(-ikr)/4\pi r \left(-(1+ikr)/r\right) \cdot dr/dn \f$ */
class helmholtz_3d_H_kernel;

/** \brief traits of the helmholtz H kernel */
template<>
struct kernel_traits<helmholtz_3d_H_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d>, normal_jacobian<space_3d> >::type trial_input_t;
	/** \brief the kernel output type */
	typedef helmholtz_3d_h_wall<space_3d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef std::complex<space_3d::scalar_t> result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = false;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 2;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
};


/** \brief 3D Helmholtz G kernel \f$ 1/4\pi r\f$ */
class helmholtz_3d_H_kernel :
	public kernel_base<helmholtz_3d_H_kernel>,
	public reciprocal_distance_kernel<helmholtz_3d_H_kernel>
{
public:
	using reciprocal_distance_kernel<helmholtz_3d_H_kernel>::estimate_complexity;

	helmholtz_3d_H_kernel(std::complex<double> wave_number) :
		m_wave_number(wave_number)
	{
	}

	std::complex<double> const &get_wave_number(void) const
	{
		return m_wave_number;
	}

private:
	std::complex<double> m_wave_number;
};

typedef helmholtz_3d_H_kernel helmholtz_3d_DLP_kernel;


#endif // HELMHOLTZ_KERNEL_HPP_INCLUDED


