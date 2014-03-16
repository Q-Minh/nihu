// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
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

class poisson_ratio_data
{
public:
	poisson_ratio_data(double nu) :	m_nu(nu) {}
	double get_poisson_ratio(void) const { return m_nu; }

private:
	double m_nu;
};

struct Ukernel
{
	typedef Eigen::Matrix<double, 3, 3> return_type;
	
	return_type operator()(
		location_input_3d const &x,
		location_input_3d const &y,
		poisson_ratio_data const &data)
	{
		auto nu = data.get_poisson_ratio();
		auto rvec = y.get_x() - x.get_x();
		auto r = rvec.norm();
		auto gradr = rvec.normalized();
		return ( (3.-4.*nu) * return_type::Identity()
				 + (gradr * gradr.transpose()) ) / (16.*M_PI*(1.-nu)*r);
	}
};

class elastostatics_3d_U_kernel;

template <>
struct kernel_traits<elastostatics_3d_U_kernel>
{
	typedef location_input_3d test_input_t;
	typedef location_input_3d trial_input_t;
	typedef collect<poisson_ratio_data> data_t;
	typedef single_brick_wall<Ukernel>::type output_t;
	typedef gauss_family_tag quadrature_family_t;
	static bool const is_symmetric = true;
	typedef asymptotic::inverse<1> far_field_behaviour_t;
	static bool const is_singular = true;
	typedef interval_estimator<
		typename reciprocal_distance_kernel_interval<1, GLOBAL_ACCURACY>::type
	> complexity_estimator_t;
};

template <>
struct singular_kernel_traits<elastostatics_3d_U_kernel>
{
	typedef asymptotic::inverse<1> singularity_type_t;
	static unsigned const singular_quadrature_order = 7;
	typedef elastostatics_3d_U_kernel singular_kernel_ancestor_t;
};

class elastostatics_3d_U_kernel :
	public kernel_base<elastostatics_3d_U_kernel>
{
public:
	elastostatics_3d_U_kernel(double nu) :
		kernel_base<elastostatics_3d_U_kernel>(poisson_ratio_data(nu)) {}
};

struct Tkernel
{
	typedef Eigen::Matrix<double, 3, 3> return_type;
	
	return_type operator()(
		location_input_3d const &x,
		location_normal_input_3d const &y,
		poisson_ratio_data const &data)
	{
		auto nu = data.get_poisson_ratio();
		auto rvec = y.get_x() - x.get_x();
		auto r = rvec.norm();
		auto gradr = rvec.normalized();
		auto const &n = y.get_unit_normal();
		auto rdny = gradr.dot(n);
		return (-rdny * ( (1.-2.*nu)*return_type::Identity() + 3.*(gradr*gradr.transpose()) )
			+ (1.-2.*nu) * (gradr*n.transpose()-n*gradr.transpose())
			) / (8.*M_PI*(1.-nu)*r*r);
	}
};

class elastostatics_3d_T_kernel;

template <>
struct kernel_traits<elastostatics_3d_T_kernel>
{
	typedef location_input_3d test_input_t;
	typedef location_normal_input_3d trial_input_t;
	typedef collect<poisson_ratio_data> data_t;
	typedef single_brick_wall<Tkernel>::type output_t;
	typedef gauss_family_tag quadrature_family_t;
	static bool const is_symmetric = false;
	typedef asymptotic::inverse<2> far_field_behaviour_t;
	static bool const is_singular = true;
	typedef interval_estimator<
		typename reciprocal_distance_kernel_interval<2, GLOBAL_ACCURACY>::type
	> complexity_estimator_t;
};

template <>
struct singular_kernel_traits<elastostatics_3d_T_kernel>
{
	typedef asymptotic::inverse<2> singularity_type_t;
	static unsigned const singular_quadrature_order = 7;
	typedef elastostatics_3d_T_kernel singular_kernel_ancestor_t;
};

class elastostatics_3d_T_kernel :
	public kernel_base<elastostatics_3d_T_kernel>
{
public:
	elastostatics_3d_T_kernel(double nu) :
		kernel_base<elastostatics_3d_T_kernel>(poisson_ratio_data(nu))
	{
	}
};


#include "guiggiani_1992.hpp"

/** \brief specialisation of class ::polar_laurent_coeffs for the ::elastostatics_3d_U_kernel */
template <>
class polar_laurent_coeffs<elastostatics_3d_U_kernel>
{
public:
	template <class guiggiani>
	static void eval(guiggiani &obj)
	{
		obj.m_Fcoeffs[0].setZero();
		obj.m_Fcoeffs[1].setZero();
	}
};

/** \brief specialisation of class ::polar_laurent_coeffs for the ::elastostatics_3d_T_kernel */
template <>
class polar_laurent_coeffs<elastostatics_3d_T_kernel>
{
public:
	template <class guiggiani>
	static void eval(guiggiani &obj)
	{
        auto const &jac0 = obj.m_Jvec_series[0];
        auto d0vec = obj.m_rvec_series[0].normalized();
        auto nu = obj.m_kernel.get_data().get_poisson_ratio();

        Eigen::Matrix<double, 3, 3> res = (d0vec*jac0.transpose()) - (jac0*d0vec.transpose());
		res *= (1.-2.*nu)/(1.-nu)/obj.m_A/obj.m_A/(8.*M_PI);

		obj.m_Fcoeffs[0] = block_product(Eigen::Matrix<double, 1, 1>(), res, obj.m_N_series[0]);
		obj.m_Fcoeffs[1].setZero();
	}
};


#endif // ELASTOSTATICS_KERNEL_HPP_INCLUDED

