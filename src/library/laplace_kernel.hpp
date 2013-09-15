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
 * \file laplace_kernel.hpp
 * \ingroup library
 * \brief implementation of kernels of the laplace equation \f$ \nabla^2 p = 0 \f$
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef LAPLACE_KERNEL_HPP_INCLUDED
#define LAPLACE_KERNEL_HPP_INCLUDED

#include <cmath>
#include "../core/global_definitions.hpp"
#include "../core/kernel.hpp"
#include "../core/gaussian_quadrature.hpp"
#include "../util/collection.hpp"
#include "interval_estimator.hpp"
#include "location_normal.hpp"
#include "basic_bricks.hpp"
#include "reciprocal_kernel_intervals.hpp"


/** \brief a brick representing a 2D laplace kernel \f$ -\log r /2\pi \f$
 * \tparam scalar scalar type of the coordinate space the distance is defined over
 */
template <class scalar>
struct laplace_2d_g_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef scalar result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel type
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
			m_laplace_g(-std::log(wall::get_distance()) / (2.0 * M_PI))
		{
		}

		/** \brief return laplace g kernel
		 * \return laplace g kernel
		 */
		scalar const & get_laplace_g(void) const
		{
			return m_laplace_g;
		}

		/** \brief return laplace g kernel
		 * \return laplace g kernel
		 */
		scalar const & get_result(void) const
		{
			return m_laplace_g;
		}

	private:
		scalar m_laplace_g;
	};
};


/** \brief combination of ::distance_vector_brick, ::distance_brick and ::laplace_2d_g_brick into a wall
 * \tparam space the coordinate space the laplace kernel is defined over
 */
template <class scalar>
struct laplace_2d_g_wall : build<
	distance_vector_brick<space<scalar, 2> >,
	distance_brick<scalar>,
	laplace_2d_g_brick<scalar>
> {};


// forward declaration
class laplace_2d_SLP_kernel;

/** \brief traits of the laplace 2D G kernel */
template<>
struct kernel_traits<laplace_2d_SLP_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_2d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_2d> >::type trial_input_t;
	/** \brief the data type */
	typedef collect<empty_data> data_t;
	/** \brief the kernel output type */
	typedef laplace_2d_g_wall<space_2d::scalar_t>::type output_t;
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

/** \brief 2D laplace kernel \f$ -\ln r/2\pi \f$ */
class laplace_2d_SLP_kernel :
	public kernel_base<laplace_2d_SLP_kernel>
{
};



/** \brief a brick representing a 2D laplace derivative kernel \f$ -1/2\pi r \cdot r'_{n_y}\f$
 * \tparam scalar scalar type of the coordinate space the distance is defined over
 */
template <class scalar>
struct laplace_2d_h_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef scalar result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel type
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
			m_laplace_h(-wall::get_rdny()/wall::get_distance() / (2.0 * M_PI))
		{
		}

		/** \brief return laplace h kernel
		 * \return laplace h kernel
		 */
		scalar const & get_laplace_h(void) const
		{
			return m_laplace_h;
		}

		/** \brief return laplace h kernel
		 * \return laplace h kernel
		 */
		scalar const & get_result(void) const
		{
			return m_laplace_h;
		}

	private:
		scalar m_laplace_h;
	};
};


/** \brief combination of ::distance_vector_brick, ::distance_brick, ::rdny_brick and ::laplace_2d_h_brick into a wall
 * \tparam space the coordinate space the laplace kernel is defined over
 */
template <class scalar>
struct laplace_2d_h_wall : build<
	distance_vector_brick<space<scalar, 2> >,
	distance_brick<scalar>,
	rdny_brick<scalar>,
	laplace_2d_h_brick<scalar>
> {};


// forward declaration
class laplace_2d_DLP_kernel;

/** \brief traits of the laplace 2D H kernel */
template<>
struct kernel_traits<laplace_2d_DLP_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_2d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_2d>, normal_jacobian<space_2d> >::type trial_input_t;
	/** \brief the data type */
	typedef collect<empty_data> data_t;
	/** \brief the kernel output type */
	typedef laplace_2d_h_wall<space_2d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef space_2d::scalar_t result_t;
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

/** \brief 2D laplace kernel \f$ -1/2\pi r \cdot r'_{n_y} \f$ */
class laplace_2d_DLP_kernel :
	public kernel_base<laplace_2d_DLP_kernel>
{
};


/** \brief a brick representing a 2D laplace derivative kernel \f$ -1/2\pi r \cdot r'_{n_x}\f$
 * \tparam scalar scalar type of the coordinate space the distance is defined over
 */
template <class scalar>
struct laplace_2d_ht_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef scalar result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel type
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
			m_laplace_ht(-wall::get_rdnx()/wall::get_distance() / (2.0 * M_PI))
		{
		}

		/** \brief return laplace h kernel
		 * \return laplace ht kernel
		 */
		scalar const & get_laplace_ht(void) const
		{
			return m_laplace_ht;
		}

		/** \brief return laplace ht kernel
		 * \return laplace ht kernel
		 */
		scalar const & get_result(void) const
		{
			return m_laplace_ht;
		}

	private:
		scalar m_laplace_ht;
	};
};


/** \brief combination of ::distance_vector_brick, ::distance_brick, ::rdnx_brick and ::laplace_2d_ht_brick into a wall
 * \tparam space the coordinate space the laplace kernel is defined over
 */
template <class scalar>
struct laplace_2d_ht_wall : build<
	distance_vector_brick<space<scalar, 2> >,
	distance_brick<scalar>,
	rdnx_brick<scalar>,
	laplace_2d_ht_brick<scalar>
> {};


// forward declaration
class laplace_2d_DLPt_kernel;

/** \brief traits of the laplace 2D Ht kernel */
template<>
struct kernel_traits<laplace_2d_DLPt_kernel>
{
	/** \brief kernel trial input type */
	typedef build<location<space_2d>, normal_jacobian<space_2d> >::type test_input_t;
	/** \brief kernel test input type */
	typedef build<location<space_2d> >::type trial_input_t;
	/** \brief the data type */
	typedef collect<empty_data> data_t;
	/** \brief the kernel output type */
	typedef laplace_2d_ht_wall<space_2d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef space_2d::scalar_t result_t;
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

/** \brief 2D laplace kernel \f$ -1/2\pi r \cdot r'_{n_x} \f$ */
class laplace_2d_DLPt_kernel :
	public kernel_base<laplace_2d_DLPt_kernel>
{
};



/** \brief a brick representing a 2D laplace hypersingular kernel \f$ \dots \f$
 * \tparam scalar scalar type of the coordinate space the distance is defined over
 */
template <class scalar>
struct laplace_2d_hyper_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef scalar result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel type
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
			m_laplace_hyper(
				1.0/(2.0 * M_PI)/wall::get_distance()/wall::get_distance() * (
					test_input.get_unit_normal().dot(trial_input.get_unit_normal()) +
					2.0 * wall::get_rdny()*wall::get_rdnx()))
		{
		}

		/** \brief return laplace hypersingular kernel
		 * \return laplace hypersingular kernel
		 */
		result_t const & get_laplace_hyper(void) const
		{
			return m_laplace_hyper;
		}

		/** \brief return laplace hypersingular kernel
		 * \return laplace hypersingular kernel
		 */
		result_t const & get_result(void) const
		{
			return m_laplace_hyper;
		}

	private:
		result_t m_laplace_hyper;
	};
};


/** \brief combination of several bricks into a laplace_2d_hyper wall
 * \tparam space the coordinate space the laplace kernel is defined over
 */
template <class scalar>
struct laplace_2d_hyper_wall : build<
	distance_vector_brick<space<scalar, 2> >,
	distance_brick<scalar>,
	rdny_brick<scalar>,
	rdnx_brick<scalar>,
	laplace_2d_hyper_brick<scalar>
> {};


// forward declaration
class laplace_2d_HSP_kernel;

/** \brief traits of the laplace 2D Hypersingular kernel */
template<>
struct kernel_traits<laplace_2d_HSP_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_2d>, normal_jacobian<space_2d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_2d>, normal_jacobian<space_2d>  >::type trial_input_t;
	/** \brief the data type */
	typedef collect<empty_data> data_t;
	/** \brief the kernel output type */
	typedef laplace_2d_hyper_wall<space_2d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef typename output_t::result_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 2;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
	/** \brief the kernel complexity estimator class */
	typedef interval_estimator<
		typename reciprocal_distance_kernel_interval<singularity_order, GLOBAL_ACCURACY>::type
	> complexity_estimator_t;
};

/** \brief 2D laplace kernel \f$ \dots \f$ */
class laplace_2d_HSP_kernel :
	public kernel_base<laplace_2d_HSP_kernel>
{
};




/** \brief a brick representing a 3D laplace kernel \f$ 1/4\pi r \f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct laplace_3d_g_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef scalar result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel type
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
			m_laplace_g(1.0 / wall::get_distance() / (4.0 * M_PI))
		{
		}

		/** \brief return laplace g kernel
		 * \return laplace g kernel
		 */
		scalar const & get_laplace_g(void) const
		{
			return m_laplace_g;
		}

		/** \brief return laplace g kernel
		 * \return laplace g kernel
		 */
		scalar const & get_result(void) const
		{
			return m_laplace_g;
		}

	private:
		scalar m_laplace_g;
	};
};


/** \brief combination of ::distance_vector_brick, ::distance_brick and ::laplace_3d_g_brick into a wall
 * \tparam space the coordinate space the laplace kernel is defined over
 */
template <class scalar>
struct laplace_3d_g_wall : build<
	distance_vector_brick<space<scalar, 3> >,
	distance_brick<scalar>,
	laplace_3d_g_brick<scalar>
> {};


// forward declaration
class laplace_3d_SLP_kernel;

/** \brief traits of the laplace 3D G kernel */
template<>
struct kernel_traits<laplace_3d_SLP_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d> >::type trial_input_t;
	/** \brief the data type */
	typedef collect<empty_data> data_t;
	/** \brief the kernel output type */
	typedef laplace_3d_g_wall<space_3d::scalar_t>::type output_t;
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

/** \brief Single layer potential kernel of the laplace equation in 3D \f$ 1/4\pi r\f$ */
class laplace_3d_SLP_kernel :
	public kernel_base<laplace_3d_SLP_kernel>
{
};


/** \brief a brick representing a laplace derivative kernel \f$ -1/4\pi r^2 \cdot r'_{n_y} \f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct laplace_3d_h_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef scalar result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel type
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
			m_laplace_h(-wall::get_laplace_g() * wall::get_rdny() / wall::get_distance())
		{
		}

		/** \brief return laplace h kernel
		 * \return laplace h kernel
		 */
		scalar const &get_laplace_h(void) const
		{
			return m_laplace_h;
		}

		/** \brief return laplace h kernel
		 * \return laplace h kernel
		 */
		scalar const & get_result(void) const
		{
			return get_laplace_h();
		}

	private:
		scalar m_laplace_h;
	};
};


/** \brief combination of laplace_g_wall and laplace_h_brick into a wall
 * \tparam space the coordinate space the laplace kernel is defined over
 */
template <class scalar>
struct laplace_3d_h_wall : build<
	distance_vector_brick<space<scalar, 3> >,
	distance_brick<scalar>,
	rdny_brick<scalar>,
	laplace_3d_g_brick<scalar>,
	laplace_3d_h_brick<scalar>
> {};


// forward declaration
class laplace_3d_DLP_kernel;

/** \brief traits of the laplace H kernel */
template<>
struct kernel_traits<laplace_3d_DLP_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d>, normal_jacobian<space_3d> >::type trial_input_t;
	/** \brief the data type */
	typedef collect<empty_data> data_t;
	/** \brief the kernel output type */
	typedef laplace_3d_h_wall<space_3d::scalar_t>::type output_t;
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

/** \brief Double layer potential kernel of the Laplace equation in 3D \f$ -1/4\pi r^2 \cdot dr/dn \f$ */
class laplace_3d_DLP_kernel :
	public kernel_base<laplace_3d_DLP_kernel>
{
};


/** \brief a brick representing a Laplace derivative kernel \f$ -1/4\pi r^2 \cdot r'_{n_x} \f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct laplace_3d_ht_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef scalar result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel type
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
			m_laplace_ht(-wall::get_laplace_g() * wall::get_rdnx() / wall::get_distance())
		{
		}

		/** \brief return laplace ht kernel
		 * \return laplace ht kernel
		 */
		scalar const &get_laplace_ht(void) const
		{
			return m_laplace_ht;
		}

		/** \brief return laplace ht kernel
		 * \return laplace ht kernel
		 */
		scalar const & get_result(void) const
		{
			return get_laplace_ht();
		}

	private:
		scalar m_laplace_ht;
	};
};


/** \brief combination of ::distance_vector_brick, ::distance_brick, ::rdnx_brick, laplace_3d_g_brick and laplace_3d_ht_brick into a wall
 * \tparam scalar the scalar type of the laplace ht kernel result
 */
template <class scalar>
struct laplace_3d_ht_wall : build<
	distance_vector_brick<space<scalar, 3> >,
	distance_brick<scalar>,
	rdnx_brick<scalar>,
	laplace_3d_g_brick<scalar>,
	laplace_3d_ht_brick<scalar>
> {};


// forward declaration
class laplace_3d_DLPt_kernel;

/** \brief traits of the laplace H kernel */
template<>
struct kernel_traits<laplace_3d_DLPt_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_3d>, normal_jacobian<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d> >::type trial_input_t;
	/** \brief the data type */
	typedef collect<empty_data> data_t;
	/** \brief the kernel output type */
	typedef laplace_3d_ht_wall<space_3d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef space_3d::scalar_t result_t;
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

/** \brief 3D laplace derivative kernel \f$ -1/4\pi r^2 \cdot r'_{n_x} \f$ */
class laplace_3d_DLPt_kernel :
	public kernel_base<laplace_3d_DLPt_kernel>
{
};


/** \brief a brick representing a laplace double derivative kernel \f$ 1/4\pi r^3 \cdot \left( n_x n_y + 3 r'_{n_x} r'_{n_y} \right) \f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct laplace_3d_hyper_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef scalar result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel type
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
			m_laplace_hyper(
				wall::get_laplace_g() / wall::get_distance() / wall::get_distance() * (
					test_input.get_unit_normal().dot(trial_input.get_unit_normal()) +
					3.0 * wall::get_rdnx() * wall::get_rdny()
				))
		{
		}

		/** \brief return laplace ht kernel
		 * \return laplace ht kernel
		 */
		scalar const &get_laplace_hyper(void) const
		{
			return m_laplace_hyper;
		}

		/** \brief return laplace ht kernel
		 * \return laplace ht kernel
		 */
		scalar const & get_result(void) const
		{
			return get_laplace_hyper();
		}

	private:
		scalar m_laplace_hyper;
	};
};


/** \brief combination of ::distance_vector_brick, ::distance_brick, ::rdnx_brick, laplace_3d_g_brick and laplace_3d_ht_brick into a wall
 * \tparam scalar the scalar type of the laplace ht kernel result
 */
template <class scalar>
struct laplace_3d_hyper_wall : build<
	distance_vector_brick<space<scalar, 3> >,
	distance_brick<scalar>,
	rdnx_brick<scalar>,
	rdny_brick<scalar>,
	laplace_3d_g_brick<scalar>,
	laplace_3d_hyper_brick<scalar>
> {};


// forward declaration
class laplace_3d_HSP_kernel;

/** \brief traits of the laplace H kernel */
template<>
struct kernel_traits<laplace_3d_HSP_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_3d>, normal_jacobian<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d>, normal_jacobian<space_3d> >::type trial_input_t;
	/** \brief the data type */
	typedef collect<empty_data> data_t;
	/** \brief the kernel output type */
	typedef laplace_3d_hyper_wall<space_3d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef space_3d::scalar_t result_t;
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

/** \brief 3D laplace derivative kernel \f$ 1/4\pi r^3 \cdot \left( n_x n_y + 3 r'_{n_x} r'_{n_y} \right) \f$ */
class laplace_3d_HSP_kernel :
	public kernel_base<laplace_3d_HSP_kernel>
{
};

#endif // LAPLACE_KERNEL_HPP_INCLUDED

