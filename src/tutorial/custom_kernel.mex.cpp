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

#ifndef ELASTOSTATICS_KERNEL_HPP_INCLUDED
#define ELASTOSTATICS_KERNEL_HPP_INCLUDED

#include "../core/global_definitions.hpp"
#include "../core/kernel.hpp"
#include "../core/gaussian_quadrature.hpp"
#include "../util/collection.hpp"
#include "location_normal.hpp"
#include "basic_bricks.hpp"


//![parameter]
class poisson_ratio_data
{
public:
	poisson_ratio_data(double nu) :	m_nu(nu) {}
	double get_poisson_ratio(void) const { return m_nu; }

private:
	double m_nu;
};
//![parameter]

//![ufunctor]
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
		return ( (3.-4.*nu) * return_type::Identity() + (gradr * gradr.transpose()) ) / (16.*M_PI*(1.-nu)*r);
	}
};
//![ufunctor]

//![udeclare]
class elastostatics_3d_U_kernel;
//![udeclare]

//![utraits]
template <>
struct kernel_traits<elastostatics_3d_U_kernel>
{
	typedef location_input_3d test_input_t;
	typedef location_input_3d trial_input_t;
	typedef collect<poisson_ratio_data> data_t;
	typedef single_brick_wall<Ukernel>::type output_t;
	enum { result_dimension = 3 };
	typedef gauss_family_tag quadrature_family_t;
	static bool const is_symmetric = true;
	typedef asymptotic::inverse<1> far_field_behaviour_t;
	static bool const is_singular = true;
};
//![utraits]

//![usingtraits]
template <>
struct singular_kernel_traits<elastostatics_3d_U_kernel>
{
	typedef asymptotic::inverse<1> singularity_type_t;
	static unsigned const singular_quadrature_order = 7;
};
//![usingtraits]

//![udefine]
class elastostatics_3d_U_kernel : public kernel_base<elastostatics_3d_U_kernel>
{
public:
	elastostatics_3d_U_kernel(double nu) :
		kernel_base<elastostatics_3d_U_kernel>(poisson_ratio_data(nu))
	{
	}
};
//![udefine]

//![tkernel]
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
	enum { result_dimension = 3 };
	typedef gauss_family_tag quadrature_family_t;
	static bool const is_symmetric = false;
	typedef asymptotic::inverse<2> far_field_behaviour_t;
	static bool const is_singular = true;
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
//![tkernel]


#include "guiggiani_1992.hpp"

template <>
class polar_laurent_coeffs<elastostatics_3d_T_kernel>
{
public:
	template <class guiggiani>
	static void eval(guiggiani &obj)
	{
        auto const &jac0 = obj.template get_Jvec_series<0>();
        auto d0vec = obj.template get_rvec_series<0>().normalized();
        auto nu = obj.get_kernel_data().get_poisson_ratio();

        Eigen::Matrix<double, 3, 3> res = (d0vec*jac0.transpose()) - (jac0*d0vec.transpose());
		res *= (1.-2.*nu)/(1.-nu)/(8.*M_PI);

		obj.template set_laurent_coeff<0>(semi_block_product(res, obj.template get_shape_series<0>()));
	}
};


#endif // ELASTOSTATICS_KERNEL_HPP_INCLUDED

