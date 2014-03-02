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
 * \file elastostatics_kernel.hpp
 * \ingroup library
 * \brief implementation of fundamental solutions of elastostatics
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef ELASTOSTATICS_KERNEL_HPP_INCLUDED
#define ELASTOSTATICS_KERNEL_HPP_INCLUDED

#include "../core/global_definitions.hpp"
#include "../core/kernel.hpp"
#include "../core/gaussian_quadrature.hpp"
#include "../util/collection.hpp"
#include "interval_estimator.hpp"
#include "location_normal.hpp"
#include "basic_bricks.hpp"
#include "reciprocal_kernel_intervals.hpp"

/** \brief kernel data that stores the Poisson's ratio */
template <class poisson_type>
class poisson_ratio_data
{
public:
	/** \brief constructor setting the Poisson's ratio */
	poisson_ratio_data(poisson_type const &nu = poisson_type()) :
		m_poisson_ratio(nu)
	{
	}

	/** \brief return Poisson's ratio */
	poisson_type const &get_poisson_ratio(void) const { return m_poisson_ratio; }

private:
	poisson_type m_poisson_ratio;
};

/** \brief a brick representing a 3D displacement kernel of elastostatics
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct elastostatics_3d_U_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef Eigen::Matrix<scalar, 3, 3> result_t;

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
			wall(test_input, trial_input, kernel_data)
		{
			auto const &nu = kernel_data.get_poisson_ratio();
			auto const &rvec = wall::get_distance_vector();
			auto const &r = wall::get_distance();
			auto gradr = rvec / r;
			m_U = (3.0 - 4.0*nu) * result_t::Identity() + (gradr * gradr.transpose());
			m_U /= (16.0*M_PI*(1.0-nu)*r);
		}

		/** \brief return U kernel */
		result_t const & get_U(void) const { return m_U; }

		/** \brief return U kernel */
		result_t const & get_result(void) const { return m_U; }

	private:
		result_t m_U;
	};
};

/** \brief combination of common bricks and elastostatics_3d_U_brick into a wall
 * \tparam scalar the scalar type of the coordinate space
 */
template <class scalar>
struct elastostatics_3d_U_wall : build<
	distance_vector_brick<space<scalar, 3> >,
	distance_brick<scalar>,
	elastostatics_3d_U_brick<scalar>
> {};


// forward declaration
class elastostatics_3d_U_kernel;

/** \brief traits of the 3D U kernel in elastostatics */
template <>
struct kernel_traits<elastostatics_3d_U_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d> >::type trial_input_t;
	/** \brief kernel data type */
	typedef collect<poisson_ratio_data<double> > data_t;
	/** \brief the kernel output type */
	typedef elastostatics_3d_U_wall<space_3d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef output_t::result_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief indicates if kernel is singular */
	static bool const is_singular = true;

	/** \brief the far field asymptotic behaviour of the kernel */
	typedef asymptotic::inverse<1> far_field_behaviour_t;
	/** \brief the kernel complexity estimator class */
	typedef interval_estimator<
		typename reciprocal_distance_kernel_interval<1, GLOBAL_ACCURACY>::type
	> complexity_estimator_t;
};

/** \brief singular traits of the 3D U kernel in elastostatics */
template <>
struct singular_kernel_traits<elastostatics_3d_U_kernel>
{
	/** \brief kernel singularity type */
	typedef asymptotic::inverse<1> singularity_type_t;
	/** \brief quadrature order used to generate blind singular quadratures */
	static unsigned const singular_quadrature_order = 7;
};


/** \brief 3D displacement kernel in elastostatics */
class elastostatics_3d_U_kernel :
	public kernel_base<elastostatics_3d_U_kernel>
{
public:
	elastostatics_3d_U_kernel(double nu) :
		kernel_base<elastostatics_3d_U_kernel>(poisson_ratio_data<double>(nu))
	{
	}
};

#endif // ELASTOSTATICS_KERNEL_HPP_INCLUDED

