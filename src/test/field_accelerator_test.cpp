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

#include "../bem/field_type_accelerator.hpp"
#include "../bem/gaussian_quadrature.hpp"

typedef quad_1_elem elem_t;
typedef field<elem_t, isoparametric_field, function_field> field_t;

#include <iostream>

int main(void)
{
	typedef gauss_family_tag gauss;
	
	field_type_accelerator<field_t, gauss> fta(5);
	for (auto it = fta.cbegin(); it != fta.cend(); ++it)
	{
		std::cout << it->get_quadrature_elem().get_xi().transpose() << std::endl;
		std::cout << it->get_shape().transpose() << std::endl;
	}
	return 0;
}
