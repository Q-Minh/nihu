#include "core/mesh.hpp"
#include "core/element.hpp"
#include "core/element_match.hpp"
#include "core/function_space.hpp"

#include <iostream>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

typedef tmp::vector<tria_2_elem, quad_2_elem> elem_type_vector_t;
typedef mesh<elem_type_vector_t> mesh_t;

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
		quad_2_elem::id, 0, 1, 2,  5,  8, 7, 6, 3, 4,
		tria_2_elem::id, 11, 9, 6, 7, 8, 10, 0, 0, 0,
		tria_2_elem::id, 8, 12, 14, 13, 11, 10, 0, 0, 0;
		
	// Create a mesh
	mesh_t msh(nodes, elements);
	
	auto const & test_space = isoparametric_view(msh);
	
	typedef function_space_view<mesh_t, field_option::isoparametric>::field_type_vector_t field_type_vector_t;

	typedef tmp::deref<tmp::begin<field_type_vector_t>::type >::type field_1_t;
	typedef tmp::at<field_type_vector_t, tmp::integer<int, 1> >::type field_2_t;
	
	auto begin_1 = test_space.field_begin<field_1_t>();
//	auto begin_2 = test_space.field_begin<field_2_t>();
	
	auto it = begin_1; ++it;
	auto match = element_match_eval(*begin_1, *(it));
	
	auto overlap = match.get_overlap();
	
	std::cout << "N: " << overlap.get_num() << std::endl;
	std::cout << "1: " << overlap.get_ind1() << std::endl;
	std::cout << "2: " << overlap.get_ind2() << std::endl;
	
	return 0;
}