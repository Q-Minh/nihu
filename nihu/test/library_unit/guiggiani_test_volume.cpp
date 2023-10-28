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

#include <iostream>

#include "library/guiggiani_1992.hpp"
#include "library/lib_element.hpp"
#include "library/location_normal.hpp"
#include "library/laplace_kernel.hpp"

typedef NiHu::quad_1_volume_elem Elem;
typedef NiHu::field_view<Elem, NiHu::field_option::constant> Field;


class MyKernel;

namespace NiHu
{
template <>
struct kernel_traits<MyKernel>
{
	typedef location_input_2d test_input_t;
	typedef location_input_2d trial_input_t;
	typedef double result_t;
	enum { result_rows = 1, result_cols = 1 };
	typedef gauss_family_tag quadrature_family_t;
	static bool const is_symmetric = true;
	typedef asymptotic::inverse<2> far_field_behaviour_t;
	static bool const is_singular = true;
};

/** \brief the singular properties of the kernel */
template <>
struct singular_kernel_traits<MyKernel>
{
	typedef asymptotic::inverse<2> singularity_type_t;
	static unsigned const singular_quadrature_order = 7;
	typedef MyKernel singular_core_t;
};
}

namespace NiHu
{
template <>
class polar_laurent_coeffs<MyKernel>
{
public:
	template <class guiggiani>
	static void eval(guiggiani &obj)
	{
        auto const &J0 = obj.get_J_series(_0());
        auto const &N0 = obj.get_shape_series(_0());

		double res = 1.;

		obj.set_laurent_coeff(_m1(), semi_block_product(res, N0*J0));
	}
};
}

class MyKernel : public NiHu::kernel_base<MyKernel>
{
public:
	double operator()(
		NiHu::location_input_2d const &x, NiHu::location_input_2d const &y) const
	{
		return 1./(x.get_x()-y.get_x()).squaredNorm();
	}
};

int main(void)
{
	Elem::coords_t coords;
	coords << 0., 1., 1., 0.,
	 0., 0., 1., 1.;
	Elem elem(coords);

	MyKernel kernel;

	NiHu::guiggiani<Field, MyKernel, 5> gui(elem, kernel);
	Eigen::Matrix<double, 1, 1> Res;
	gui.integrate(Res, Elem::domain_t::get_center());

	return 0;
}

