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
#include "../bem/mesh.hpp"
#include "../bem/singularity_check.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_vector;
typedef Mesh<elem_vector> mesh_t;

typedef quad_1_elem elem_t;
typedef field<elem_t, constant_field, dirac_field> test_field_t;
typedef field<elem_t, constant_field, function_field> trial_field_t;

class kernel {};

int main(void)
{

	double c[] = {
		0, 0, 0,
		1, 0, 0,
		2, 0, 0,
		0, 1, 0,
		1, 1, 0,
		2, 1, 0,
		3, 1, 0
	};
	unsigned e[] = {
		4,  0, 3, 4, 1,
		4,  5, 2, 1, 4,
		3,  2, 6, 5, 0
	};

	mesh_t mesh;
	for (unsigned i = 0; i < 7; ++i)
		mesh.add_node(c+i*3);
	for (unsigned i = 0; i < 3; ++i)
		mesh.add_elem(e+i*5);

	tria_1_elem e1 = *(mesh.begin<tria_1_elem>());
	quad_1_elem e2 = *(mesh.begin<quad_1_elem>());
	quad_1_elem e3 = *(mesh.begin<quad_1_elem>()+1);

	elem_t test_elem = elem_t(elem_t::coords_t());
	test_field_t test_field(test_elem);
	trial_field_t trial_field(test_elem);

	std::cout <<
		(singularity_check<kernel, test_field_t, trial_field_t>::eval(test_field, trial_field) == FACE_MATCH)
		<< std::endl;

	element_overlapping o1 = e1.get_overlapping(e1);
	element_overlapping o2 = e2.get_overlapping(e1);
	element_overlapping o3 = e2.get_overlapping(e3);
	element_overlapping o4 = e1.get_overlapping(e2);

	std::cout << "Overlapping tria  (2,6,5)   -> quad1 (0,1,4,3) : " <<
		o4.get_num() << " " << o4.get_ind1() << " " << o4.get_ind2() << std::endl;
	std::cout << "Overlapping tria  (2,6,5)   -> tria  (2,6,5)   : " <<
		o1.get_num() << " " << o1.get_ind1() << " " << o1.get_ind2() << std::endl;
	std::cout << "Overlapping quad1 (0,1,4,3) -> tria  (2,6,5)   : " <<
		o2.get_num() << " " << o2.get_ind1() << " " << o2.get_ind2() << std::endl;
	std::cout << "Overlapping quad1 (0,3,4,1) -> quad2 (5,2,1,4) : " <<
		o3.get_num() << " " << o3.get_ind1() << " " << o3.get_ind2() << std::endl;

	return 0;
}
