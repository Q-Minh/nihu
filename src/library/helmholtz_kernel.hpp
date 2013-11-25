// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2013  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2013  Peter Rucz <rucz@hit.bme.hu>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

/**
 * \file helmholtz_kernel.hpp
 * \ingroup library
 * \brief implementation of kernels of the Helmholtz equation \f$ \nabla^2 p + k^2 p = 0 \f$
 */

#ifndef HELMHOLTZ_KERNEL_HPP_INCLUDED
#define HELMHOLTZ_KERNEL_HPP_INCLUDED

#include <cmath>
#include <complex>

#include "../core/global_definitions.hpp"

#include "../core/kernel.hpp"
#include "../core/gaussian_quadrature.hpp"

#include "location_normal.hpp"

#include "basic_bricks.hpp"
#include "../util/collection.hpp"
#include "../util/math_functions.hpp"

#include "reciprocal_kernel_intervals.hpp"
#include "interval_estimator.hpp"

/**
 * \brief kernel data that stores the wave number
 * \tparam wave_number_type the wave number type
 */
template <class wave_number_type>
class wave_number_data
{
public:
	/** \brief constructor setting the wave number
	 * \param [in] wn the wave number to set
	 */
	wave_number_data(wave_number_type const &wn = wave_number_type()) :
		m_wave_number(wn)
	{
	}

	/** \brief return wave number
	 * \return wave number
	 */
	wave_number_type const &get_wave_number(void) const
	{
		return m_wave_number;
	}

private:
	wave_number_type m_wave_number;
};

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
		/** \brief the result type */
		typedef std::complex<scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel data type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel_data the kernel data instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_data_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_data_t const &kernel_data) :
			wall(test_input, trial_input, kernel_data),
			m_ikr(std::complex<scalar>(0.0,1.0) * kernel_data.get_wave_number() * wall::get_distance())
		{
		}

		/** \brief return ikr
		 * \return ikr
		 */
		result_t const &get_ikr(void) const
		{
			return m_ikr;
		}

	private:
		result_t m_ikr;
	};
};


/** \brief a brick representing a 2D Helmholtz kernel \f$ -i/4 H_0(kr) \f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct helmholtz_2d_g_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef std::complex<scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel data type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel_data the kernel data instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_data_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_data_t const &kernel_data) :
			wall(test_input, trial_input, kernel_data),
			m_helmholtz_g(std::complex<scalar>(0.0, -.25) *
				bessel::H<0>(kernel_data.get_wave_number()*wall::get_distance()))
		{
		}

		/** \brief return Helmholtz g kernel
		 * \return Helmholtz g kernel
		 */
		result_t const &get_helmholtz_g(void) const
		{
			return m_helmholtz_g;
		}

		/** \brief return Helmholtz g kernel
		 * \return Helmholtz g kernel
		 */
		result_t const &get_result(void) const
		{
			return m_helmholtz_g;
		}

	private:
		result_t m_helmholtz_g;
	};
};


/** \brief a brick representing a 3D Helmholtz kernel \f$ exp(-ikr) / 4\pi r \f$
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
		/** \brief the result type */
		typedef std::complex<scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel data type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel_data the kernel data instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_data_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_data_t const &kernel_data) :
			wall(test_input, trial_input, kernel_data),
			m_helmholtz_g(std::exp(-wall::get_ikr()) / wall::get_distance() / (4.0 * M_PI))
		{
		}

		/** \brief return Helmholtz g kernel
		 * \return Helmholtz g kernel
		 */
		result_t const &get_helmholtz_g(void) const
		{
			return m_helmholtz_g;
		}

		/** \brief return Helmholtz g kernel
		 * \return Helmholtz g kernel
		 */
		result_t const &get_result(void) const
		{
			return m_helmholtz_g;
		}

	private:
		result_t m_helmholtz_g;
	};
};


/** \brief combination of ::distance_vector_brick, ::distance_brick and ::helmholtz_2d_g_brick into a wall
 * \tparam space the coordinate space the Helmholtz kernel is defined over
 */
template <class scalar>
struct helmholtz_2d_g_wall : build<
	distance_vector_brick<space<scalar, 2> >,
	distance_brick<scalar>,
	helmholtz_2d_g_brick<scalar>
> {};

/** \brief combination of ::distance_vector_brick, ::distance_brick, ::ikr_brick and ::helmholtz_3d_g_brick into a wall
 * \tparam space the coordinate space the Helmholtz kernel is defined over
 */
template <class scalar>
struct helmholtz_3d_g_wall : build<
	distance_vector_brick<space<scalar, 3> >,
	distance_brick<scalar>,
	ikr_brick<scalar>,
	helmholtz_3d_g_brick<scalar>
> {};

// forward declaration
template <class wave_number_t>
class helmholtz_2d_SLP_kernel;

/** \brief traits of the Helmholtz G kernel */
template <class wave_number_t>
struct kernel_traits<helmholtz_2d_SLP_kernel<wave_number_t> >
{
	/** \brief kernel test input type */
	typedef build<location<space_2d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_2d> >::type trial_input_t;
	/** \brief kernel data type */
	typedef collect<wave_number_data<wave_number_t> > data_t;
	/** \brief the kernel output type */
	typedef helmholtz_2d_g_wall<space_2d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef typename output_t::result_t result_t;
	/** \brief the quadrature family the kernel is integrated with
	 * \todo update this quantity
	 */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) )
	 * \todo update this quantity
	 */
	static unsigned const singularity_order = 1;
	/** \brief quadrature order used to generate Duffy singular quadratures
	 * \todo update this quantity
	 */
	static unsigned const singular_quadrature_order = 7;
	/** \brief the kernel complexity estimator class
	 * \todo update this quantity
	 */
	typedef interval_estimator<
		typename reciprocal_distance_kernel_interval<singularity_order, GLOBAL_ACCURACY>::type
	> complexity_estimator_t;
};


/** \brief Single layer potential kernel of the Helmholtz equation in 2D \f$ -i/4 H_0(kr) \f$ */
template <class wave_number_t>
class helmholtz_2d_SLP_kernel :
	public kernel_base<helmholtz_2d_SLP_kernel<wave_number_t> >
{
public:
	/** \brief constructor
	 * \param [in] wave_number the wave number
	 */
	helmholtz_2d_SLP_kernel(wave_number_t const &wave_number) :
		kernel_base<helmholtz_2d_SLP_kernel<wave_number_t> >(wave_number_data<wave_number_t>(wave_number))
	{
	}
};


// forward declaration
template <class wave_number_t>
class helmholtz_3d_SLP_kernel;

/** \brief traits of the Helmholtz G kernel */
template <class wave_number_t>
struct kernel_traits<helmholtz_3d_SLP_kernel<wave_number_t> >
{
	/** \brief kernel test input type */
	typedef build<location<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d> >::type trial_input_t;
	/** \brief kernel data type */
	typedef collect<wave_number_data<wave_number_t> > data_t;
	/** \brief the kernel output type */
	typedef helmholtz_3d_g_wall<space_3d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef typename output_t::result_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 1;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
	/** \brief the kernel complexity estimator class */
	typedef interval_estimator<
		typename reciprocal_distance_kernel_interval<singularity_order, GLOBAL_ACCURACY>::type
	> complexity_estimator_t;
};


/** \brief Single layer potential kernel of the Helmholtz equation in 3D \f$ \exp(-ikr)/4\pi r\f$ */
template <class wave_number_t>
class helmholtz_3d_SLP_kernel :
	public kernel_base<helmholtz_3d_SLP_kernel<wave_number_t> >
{
public:
	/** \brief constructor
	 * \param [in] wave_number the wave number
	 */
	helmholtz_3d_SLP_kernel(wave_number_t const &wave_number) :
		kernel_base<helmholtz_3d_SLP_kernel<wave_number_t> >(wave_number_data<wave_number_t>(wave_number))
	{
	}
};


/** \brief a brick representing a 2D Helmholtz derivative kernel \f$ ik/4 H_1(kr) \cdot r'_{n_y} \f$
* \tparam scalar the scalar of the coordinate space the distance is defined over
*/
template <class scalar>
struct helmholtz_2d_h_brick
{
	/** \brief the brick template
	* \tparam the wall the brick is placed on
	*/
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef std::complex<scalar> result_t;

		/** \brief templated constructor
		* \tparam test_input_t the test input type
		* \tparam trial_input_t the trial input type
		* \tparam kernel_data_t the kernel data type
		* \param [in] test_input the test input
		* \param [in] trial_input the trial input
		* \param [in] kernel_data the kernel data instance
		*/
		template <class test_input_t, class trial_input_t, class kernel_data_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_data_t const &kernel_data) :
			wall(test_input, trial_input, kernel_data),
			m_helmholtz_h(std::complex<scalar>(0.0, .25) * kernel_data.get_wave_number() *
			bessel::H<1>(kernel_data.get_wave_number()*wall::get_distance())  * wall::get_rdny())
		{
		}

		/** \brief return Helmholtz h kernel
		* \return Helmholtz h kernel
		*/
		result_t const &get_helmholtz_h(void) const
		{
			return m_helmholtz_h;
		}

		/** \brief return Helmholtz h kernel
		* \return Helmholtz h kernel
		*/
		result_t const &get_result(void) const
		{
			return m_helmholtz_h;
		}

	private:
		result_t m_helmholtz_h;
	};
};


/** \brief a brick representing a 3D Helmholtz H kernel \f$ \exp(-ikr)/4\pi r \left(-(1+ikr)/r\right) \cdot r'_{n_y} \f$
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
		/** \brief the result type */
		typedef std::complex<scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel data type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel_data the kernel data instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_data_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_data_t const &kernel_data) :
			wall(test_input, trial_input, kernel_data),
			m_helmholtz_h(-(1.0+wall::get_ikr()) * wall::get_helmholtz_g() / wall::get_distance() * wall::get_rdny())
		{
		}

		/** \brief return helmholtz h kernel
		 * \return helmholtz h kernel
		 */
		result_t const &get_helmholtz_h(void) const
		{
			return m_helmholtz_h;
		}

		/** \brief return Helmholtz h kernel
		 * \return Helmholtz h kernel
		 */
		result_t const &get_result(void) const
		{
			return m_helmholtz_h;
		}

	private:
		result_t m_helmholtz_h;
	};
};


template <class scalar>
struct helmholtz_2d_h_wall : build<
	distance_vector_brick<space<scalar, 2> >,
	distance_brick<scalar>,
	rdny_brick<scalar>,
	helmholtz_2d_h_brick<scalar>
> {};


template <class scalar>
struct helmholtz_3d_h_wall : build<
	distance_vector_brick<space<scalar, 3> >,
	distance_brick<scalar>,
	ikr_brick<scalar>,
	helmholtz_3d_g_brick<scalar>,
	rdny_brick<scalar>,
	helmholtz_3d_h_brick<scalar>
> {};


// forward declaration
template <class wave_number_t>
class helmholtz_2d_DLP_kernel;

/** \brief traits of the Helmholtz H kernel */
template <class wave_number_t>
struct kernel_traits<helmholtz_2d_DLP_kernel<wave_number_t> >
{
	/** \brief kernel test input type */
	typedef build<location<space_2d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_2d>, normal_jacobian<space_2d> >::type trial_input_t;
	/** \brief kernel data type */
	typedef collect<wave_number_data<wave_number_t> > data_t;
	/** \brief the kernel output type */
	typedef helmholtz_2d_h_wall<space_2d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef typename output_t::result_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = false;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 1;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
	/** \brief the kernel complexity estimator class */
	typedef interval_estimator<
		typename reciprocal_distance_kernel_interval<singularity_order, GLOBAL_ACCURACY>::type
	> complexity_estimator_t;
};


/** \brief Double layer potential kernel of the Helmholtz equation in 3D \f$ \exp(-ikr)/4\pi r \left(-(1+ikr)/r\right) \cdot r'_{n_y} \f$ */
template <class wave_number_t>
class helmholtz_2d_DLP_kernel :
	public kernel_base<helmholtz_2d_DLP_kernel<wave_number_t> >
{
public:
	/** \brief constructor
	* \param [in] wave_number the wave number
	*/
	helmholtz_2d_DLP_kernel(wave_number_t const &wave_number) :
		kernel_base<helmholtz_2d_DLP_kernel<wave_number_t> >(wave_number_data<wave_number_t>(wave_number))
	{
	}
};


// forward declaration
template <class wave_number_t>
class helmholtz_3d_DLP_kernel;

/** \brief traits of the Helmholtz H kernel */
template <class wave_number_t>
struct kernel_traits<helmholtz_3d_DLP_kernel<wave_number_t> >
{
	/** \brief kernel test input type */
	typedef build<location<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d>, normal_jacobian<space_3d> >::type trial_input_t;
	/** \brief kernel data type */
	typedef collect<wave_number_data<wave_number_t> > data_t;
	/** \brief the kernel output type */
	typedef helmholtz_3d_h_wall<space_3d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef typename output_t::result_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = false;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 2;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
	/** \brief the kernel complexity estimator class */
	typedef interval_estimator<
		typename reciprocal_distance_kernel_interval<singularity_order, GLOBAL_ACCURACY>::type
	> complexity_estimator_t;
};


/** \brief Double layer potential kernel of the Helmholtz equation in 3D \f$ \exp(-ikr)/4\pi r \left(-(1+ikr)/r\right) \cdot r'_{n_y} \f$ */
template <class wave_number_t>
class helmholtz_3d_DLP_kernel :
	public kernel_base<helmholtz_3d_DLP_kernel<wave_number_t> >
{
public:
	/** \brief constructor
	 * \param [in] wave_number the wave number
	 */
	helmholtz_3d_DLP_kernel(wave_number_t const &wave_number) :
		kernel_base<helmholtz_3d_DLP_kernel<wave_number_t> >(wave_number_data<wave_number_t>(wave_number))
	{
	}
};


/** \brief a brick representing a 3D Helmholtz H kernel \f$ \exp(-ikr)/4\pi r \left(-(1+ikr)/r\right) \cdot dr/dn \f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct helmholtz_3d_ht_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef std::complex<scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel data type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel_data the kernel data instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_data_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_data_t const &kernel_data) :
			wall(test_input, trial_input, kernel_data),
			m_helmholtz_ht(-(1.0+wall::get_ikr()) * wall::get_helmholtz_g() / wall::get_distance() * wall::get_rdnx())
		{
		}

		/** \brief return Helmholtz Ht kernel
		 * \return Helmholtz Ht kernel
		 */
		result_t const &get_helmholtz_ht(void) const
		{
			return m_helmholtz_ht;
		}

		/** \brief return Helmholtz h kernel
		 * \return Helmholtz h kernel
		 */
		result_t const &get_result(void) const
		{
			return m_helmholtz_ht;
		}

	private:
		std::complex<scalar> m_helmholtz_ht;
	};
};


template <class scalar>
struct helmholtz_3d_ht_wall : build<
	distance_vector_brick<space<scalar, 3> >,
	distance_brick<scalar>,
	ikr_brick<scalar>,
	helmholtz_3d_g_brick<scalar>,
	rdnx_brick<scalar>,
	helmholtz_3d_ht_brick<scalar>
> {};


/** \brief Transposed Double layer potential kernel of the Helmholtz equation in 3D \f$ \exp(-ikr)/4\pi r \left(-(1+ikr)/r\right) \cdot r'_{n_x} \f$ */
template <class wave_number_t>
class helmholtz_3d_DLPt_kernel;

/** \brief traits of the Helmholtz Ht kernel */
template <class wave_number_t>
struct kernel_traits<helmholtz_3d_DLPt_kernel<wave_number_t> >
{
	/** \brief kernel test input type */
	typedef build<location<space_3d>, normal_jacobian<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d> >::type trial_input_t;
	/** \brief kernel data type */
	typedef collect<wave_number_data<wave_number_t> > data_t;
	/** \brief the kernel output type */
	typedef helmholtz_3d_ht_wall<space_3d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef typename output_t::result_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = false;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 2;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
	/** \brief the kernel complexity estimator class */
	typedef interval_estimator<
		typename reciprocal_distance_kernel_interval<singularity_order, GLOBAL_ACCURACY>::type
	> complexity_estimator_t;
};


/** \brief 3D Helmholtz Ht kernel \f$ \exp(-ikr)/4\pi r \left(-(1+ikr)/r\right) \cdot r'_{n_x} \f$ */
template <class wave_number_t>
class helmholtz_3d_DLPt_kernel :
	public kernel_base<helmholtz_3d_DLPt_kernel<wave_number_t> >
{
public:
	/** \brief constructor
	 * \param [in] wave_number the wave number
	 */
	helmholtz_3d_DLPt_kernel(wave_number_t const &wave_number) :
		kernel_base<helmholtz_3d_DLPt_kernel<wave_number_t> >(wave_number_data<wave_number_t>(wave_number))
	{
	}
};


/** \brief a brick representing a 3D Helmholtz H kernel \f$ \exp(-ikr)/4\pi r \left(-(1+ikr)/r\right) \cdot dr/dn \f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct helmholtz_3d_hyper_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef std::complex<scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel data type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel_data the kernel data instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_data_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_data_t const &kernel_data) :
			wall(test_input, trial_input, kernel_data),
			m_helmholtz_hyper(
				wall::get_helmholtz_g()/wall::get_distance()/wall::get_distance() * (
					(1.0 + wall::get_ikr())*test_input.get_unit_normal().dot(trial_input.get_unit_normal()) +
					(3.0 + 3.0*wall::get_ikr() + wall::get_ikr()*wall::get_ikr())*wall::get_rdny()*wall::get_rdnx()
				)
			)
		{
		}

		/** \brief return helmholtz hypersingular kernel
		 * \return helmholtz hypersingular kernel
		 */
		result_t const &get_helmholtz_hyper(void) const
		{
			return m_helmholtz_hyper;
		}

		/** \brief return Helmholtz hypersingular kernel
		 * \return Helmholtz hypersingular kernel
		 */
		result_t const &get_result(void) const
		{
			return m_helmholtz_hyper;
		}

	private:
		result_t m_helmholtz_hyper;
	};
};


template <class scalar>
struct helmholtz_3d_hyper_wall : build<
	distance_vector_brick<space<scalar, 3> >,
	distance_brick<scalar>,
	ikr_brick<scalar>,
	helmholtz_3d_g_brick<scalar>,
	rdnx_brick<scalar>,
	rdny_brick<scalar>,
	helmholtz_3d_hyper_brick<scalar>
> {};


// forward declaration
template <class wave_number_t>
class helmholtz_3d_HSP_kernel;

/** \brief traits of the Helmholtz Hyper kernel */
template <class wave_number_t>
struct kernel_traits<helmholtz_3d_HSP_kernel<wave_number_t> >
{
	/** \brief kernel test input type */
	typedef build<location<space_3d>, normal_jacobian<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d>, normal_jacobian<space_3d> >::type trial_input_t;
	/** \brief kernel data type */
	typedef collect<wave_number_data<wave_number_t> > data_t;
	/** \brief the kernel output type */
	typedef helmholtz_3d_hyper_wall<space_3d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef typename output_t::result_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 3;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
	/** \brief the kernel complexity estimator class */
	typedef interval_estimator<
		typename reciprocal_distance_kernel_interval<singularity_order, GLOBAL_ACCURACY>::type
	> complexity_estimator_t;
};


/** \brief Hypersingular kernel of the Helmholtz equation in 3D */
template <class wave_number_t>
class helmholtz_3d_HSP_kernel :
	public kernel_base<helmholtz_3d_HSP_kernel<wave_number_t> >
{
public:
	/** \brief constructor
	 * \param [in] wave_number the wave number
	 */
	helmholtz_3d_HSP_kernel(wave_number_t const &wave_number) :
		kernel_base<helmholtz_3d_HSP_kernel<wave_number_t> >(wave_number_data<wave_number_t>(wave_number))
	{
	}
};

#endif // HELMHOLTZ_KERNEL_HPP_INCLUDED

