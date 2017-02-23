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

#include "core/blind_singular_quadrature.hpp"
#include "core/singular_accelerator.hpp"
#include "library/laplace_kernel.hpp"
#include "library/lib_element.hpp"

#include <iostream>

typedef NiHu::select_singular_accelerator<
	NiHu::laplace_3d_HSP_kernel,
	NiHu::field<NiHu::quad_1_elem, NiHu::quad_1_shape_set>,
	NiHu::field<NiHu::quad_1_elem, NiHu::quad_1_shape_set>
>::type acc_t;

int main(void)
{
	std::cout << std::is_same<acc_t, NiHu::invalid_singular_accelerator>::value << std::endl;

	std::cout << std::is_same<
		NiHu::blind_transform::duffy,
		NiHu::blind_transform_selector<
            NiHu::asymptotic::inverse<1>,
            NiHu::tria_domain
		>::type
	>::value;

	return 0;
}
