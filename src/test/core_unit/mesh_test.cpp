#include <iostream>
#include "bem/mesh.hpp"

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
			
			unsigned n = 0;
			while (b!=e)
			{
				++b;
				++n;
			}
			std::cout << "# Elem (ID = " << ElemType::id << "): " << n << std::endl;
			
		}
	};
	
};


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

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
		quad_1_elem::id, 1, 2, 5, 4,
		quad_1_elem::id, 3, 4, 7, 6,
		tria_1_elem::id, 4, 5, 8, 0,
		tria_1_elem::id, 4, 8, 7, 0;

	// Create a mesh
	mesh_t msh(nodes, elements);
	
	// Call tester
	tmp::call_each<
		elem_type_vector_t,
		tester<tmp::_1>,
		const mesh_t&
	>(msh);
	
	return 0;
}