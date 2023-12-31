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

#include "element_test.h"
#include "nihu/core/element.hpp"
#include "nihu/tmp/sequence.hpp"
#include "nihu/tmp/control.hpp"
#include "nihu/tmp/vector.hpp"

#include "nihu/library/lib_element.hpp"

#include <iostream>

template <class ElemType>
struct tester
{
	struct type
	{
		void operator() (void)
		{
			std::cout << std::endl;
			std::cout << "Elem type ID: " << ElemType::id << std::endl;
			std::cout << "===================" << std::endl;

			typename ElemType::xi_t xi = ElemType::domain_t::get_center();
			std::cout << "xi center: " << xi.transpose() << std::endl;

			ElemType e(ElemType::coords_t::Random());

			typename ElemType::x_t c = e.get_center();
			std::cout << "x  center: " << c.transpose() << std::endl;
			std::cout << "nodal coords: " << std::endl << e.get_coords() << std::endl;

			typename ElemType::x_return_type x = e.get_x(xi);
			std::cout << "x(xi): " << x.transpose() << std::endl;
			typename ElemType::dx_return_type dx = e.get_dx(xi);
			std::cout << "dx(xi): " << dx.transpose() << std::endl;
			typename ElemType::ddx_return_type ddx = e.get_ddx(xi);
			std::cout << "ddx(xi): " << ddx.transpose() << std::endl;

			typename ElemType::x_t n = e.get_normal(xi);
			std::cout << "n(xi): " << n.transpose() << std::endl;

			std::cout << "Linear size estimate:\t" << e.get_linear_size_estimate() << std::endl;
		}
	};
};

int main(void)
{
	tmp::call_each<
		tmp::vector<
			NiHu::line_1_elem, NiHu::line_2_elem,
			NiHu::tria_1_elem, NiHu::quad_1_elem,
			NiHu::tria_2_elem, NiHu::quad_2_elem, NiHu::quad_28_elem
		>,
		tester<tmp::_1>
	>();

	return 0;
}

