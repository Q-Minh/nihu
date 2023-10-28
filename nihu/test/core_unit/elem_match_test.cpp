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

#include "nihu/core/mesh.hpp"
#include "nihu/core/element.hpp"
#include "nihu/core/element_match.hpp"
#include "nihu/core/function_space.hpp"
#include "nihu/library/lib_element.hpp"

#include <iostream>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

typedef tmp::vector<NiHu::tria_2_elem, NiHu::quad_2_elem> elem_type_vector_t;
typedef NiHu::mesh<elem_type_vector_t> mesh_t;

int main()
{
	// Definition of nodes
	dMatrix nodes(15,3);
	nodes <<
		 0.0,  0.0, 0.0,
		 1.0,  0.0, 0.0,
		 2.0,  0.0, 0.0,

		 0.0,  1.0, 0.0,
		 1.0,  1.0, 0.0,
		 2.0,  1.0, 0.0,

		 0.0,  2.0, 0.0,
		 1.0,  2.0, 0.0,
		 2.0,  2.0, 0.0,

		 0.0,  3.0, 0.0,
		 1.0,  3.0, 0.0,
		 0.0,  4.0, 0.0,

		 2.0,  3.0, 0.0,
		 1.0,  4.0, 0.0,
		 2.0,  4.0, 0.0;

	// Definition of elements
	uMatrix elements(3, 1+9);
	elements <<
		NiHu::quad_2_elem::id, 0, 1, 2,  5,  8, 7, 6, 3, 4,
		NiHu::tria_2_elem::id, 11, 9, 6, 7, 8, 10, 0, 0, 0,
		NiHu::tria_2_elem::id, 8, 12, 14, 13, 11, 10, 0, 0, 0;

	// Create a mesh
	mesh_t msh(nodes, elements);

	auto const & test_space = NiHu::isoparametric_view(msh);

	typedef NiHu::function_space_view<mesh_t, NiHu::field_option::isoparametric>::field_type_vector_t field_type_vector_t;

	typedef tmp::deref<tmp::begin<field_type_vector_t>::type >::type field_1_t;
//	typedef tmp::at<field_type_vector_t, tmp::integer<int, 1> >::type field_2_t;

	auto begin_1 = test_space.field_begin<field_1_t>();
//	auto begin_2 = test_space.field_begin<field_2_t>();

	auto it = begin_1; ++it;
	auto match = NiHu::element_match_eval(*begin_1, *(it));

	auto overlap = match.get_overlap();

	std::cout << "N: " << overlap.get_num() << std::endl;
	std::cout << "1: " << overlap.get_ind1() << std::endl;
	std::cout << "2: " << overlap.get_ind2() << std::endl;

	return 0;
}

