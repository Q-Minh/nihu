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

#include "core/mesh.hpp"

#include <iostream>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;
typedef tmp::vector<tria_1_elem, quad_1_elem, tria_2_elem, quad_2_elem> elem_type_vector_t;
typedef mesh<elem_type_vector_t> mesh_t;

template <class ElemType>
struct tester
{
	struct type
	{
		void operator() (mesh_t const& mesh)
		{
			auto b = mesh.template begin<ElemType>();
			auto e = mesh.template end<ElemType>();
			
			std::cout << "ID = " << ElemType::id << " : ";
			
			unsigned n = 0;
			while (b!=e)
			{
				std::cout << (*b).get_id() << ' ';
				++b;
				++n;
			}
			std::cout << "(" << n << ") " << std::endl;
		}
	};
	
};

int main(void)
{
	// Definition of nodes
	dMatrix nodes(9,3);
	nodes <<
		-1.0, -1.0, 0.0,
		 0.0, -1.0, 0.0,
		 1.0, -1.0, 0.0,

		-1.0,  0.0, 0.0,
		 0.0,  0.0, 0.0,
		 1.0,  0.0, 0.0,

		-1.0,  1.0, 0.0,
		 0.0,  1.0, 0.0,
		 1.0,  1.0, 0.0;

	// Definition of elements
	uMatrix elements(5, 1+4);
	elements <<
		quad_1_elem::id, 0, 1, 4, 3,
		tria_1_elem::id, 1, 2, 5, 0,
		quad_1_elem::id, 3, 4, 7, 6,
		tria_1_elem::id, 4, 5, 8, 0,
		tria_1_elem::id, 4, 8, 7, 0;

	// Create a mesh
	mesh_t msh(nodes, elements);
	
	std::cout << "Listing elements by type" << std::endl;
	std::cout << "========================" << std::endl;
	
	// Call tester
	tmp::call_each<
		elem_type_vector_t,
		tester<tmp::_1>,
		const mesh_t&
	>(msh);
	
	std::cout << std::endl;
	
	return 0;
}
