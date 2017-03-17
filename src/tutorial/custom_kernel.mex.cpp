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
#include "../library/location_normal.hpp"


//![parameter]
class poisson_ratio_kernel
{
public:
	poisson_ratio_kernel(double nu)
		: m_nu(nu)
	{
	}
	
	double get_poisson_ratio(void) const
	{
		return m_nu;
	}

private:
	double m_nu;
};
//![parameter]

//![udeclare]
class elastostatics_3d_U_kernel;
//![udeclare]

//![utraits]
namespace NiHu
{
template <>
struct kernel_traits<elastostatics_3d_U_kernel>
{
	typedef location_input_3d test_input_t;
	typedef location_input_3d trial_input_t;
	typedef Eigen::Matrix<double, 3, 3> result_t;
	enum { result_rows = 3, result_cols = 3 };
	typedef gauss_family_tag quadrature_family_t;
	static bool const is_symmetric = true;
	typedef asymptotic::inverse<1> far_field_behaviour_t;
	static bool const is_singular = true;
};
}
//![utraits]

//![usingtraits]
namespace NiHu
{
template <>
struct singular_kernel_traits<elastostatics_3d_U_kernel>
{
	typedef asymptotic::inverse<1> singularity_type_t;
	static unsigned const singular_quadrature_order = 7;
};
}
//![usingtraits]

//![udefine]
class elastostatics_3d_U_kernel
	: public NiHu::kernel_base<elastostatics_3d_U_kernel>
	, public poisson_ratio_kernel
{
public:
	elastostatics_3d_U_kernel(double nu) :
		poisson_ratio_kernel(nu)
	{
	}
	
	result_t operator()(
		NiHu::location_input_3d const &x,
		NiHu::location_input_3d const &y) const
	{
		auto nu = get_poisson_ratio();
		auto rvec = y.get_x() - x.get_x();
		auto r = rvec.norm();
		auto gradr = rvec.normalized();
		return ( (3.-4.*nu) * result_t::Identity() + (gradr * gradr.transpose()) ) / (16.*M_PI*(1.-nu)*r);
	}
};
//![udefine]

//![tkernel]
class elastostatics_3d_T_kernel;

namespace NiHu
{
template <>
struct kernel_traits<elastostatics_3d_T_kernel>
{
	typedef location_input_3d test_input_t;
	typedef location_normal_input_3d trial_input_t;
	typedef Eigen::Matrix<double, 3, 3> result_t;
	enum { result_rows = 3, result_cols = 3 };
	typedef gauss_family_tag quadrature_family_t;
	static bool const is_symmetric = false;
	typedef asymptotic::inverse<2> far_field_behaviour_t;
	static bool const is_singular = true;
};
}

namespace NiHu
{
template <>
struct singular_kernel_traits<elastostatics_3d_T_kernel>
{
	typedef asymptotic::inverse<2> singularity_type_t;
	static unsigned const singular_quadrature_order = 7;
	typedef elastostatics_3d_T_kernel singular_kernel_ancestor_t;
};
}

class elastostatics_3d_T_kernel
	: public NiHu::kernel_base<elastostatics_3d_T_kernel>
	, public poisson_ratio_kernel
{
public:
	elastostatics_3d_T_kernel(double nu)
		: poisson_ratio_kernel(nu)
	{
	}

	result_t operator()(
		NiHu::location_input_3d const &x,
		NiHu::location_normal_input_3d const &y) const
	{
		auto nu = get_poisson_ratio();
		auto rvec = y.get_x() - x.get_x();
		auto r = rvec.norm();
		auto gradr = rvec.normalized();
		auto const &n = y.get_unit_normal();
		auto rdny = gradr.dot(n);
		return (-rdny * ( (1.-2.*nu)*result_t::Identity() + 3.*(gradr*gradr.transpose()) )
			+ (1.-2.*nu) * (gradr*n.transpose()-n*gradr.transpose())
			) / (8.*M_PI*(1.-nu)*r*r);
	}
};
//![tkernel]


#include "../library/guiggiani_1992.hpp"

namespace NiHu
{
template <>
class polar_laurent_coeffs<elastostatics_3d_T_kernel>
{
public:
	template <class guiggiani>
	static void eval(guiggiani &obj)
	{
        auto const &r1 = obj.get_rvec_series(_1());
        auto const &j0 = obj.get_Jvec_series(_0());
        auto const &N0 = obj.get_shape_series(_0());
        auto nu = obj.get_poisson_ratio();
        Eigen::Matrix<double, 3, 3> res = ((r1*j0.transpose())-(j0*r1.transpose()))
			* (1.-2.*nu)/(1.-nu)/(8.*M_PI);
		obj.set_laurent_coeff(_m1(), semi_block_product(res, N0));
	}
};
}

#endif // ELASTOSTATICS_KERNEL_HPP_INCLUDED

