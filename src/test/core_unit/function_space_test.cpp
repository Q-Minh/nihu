#include "bem/function_space.hpp"


typedef tmp::vector<tria_1_elem, quad_1_elem> elem_type_vector_t;
typedef mesh<elem_type_vector_t> mesh_t;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;


template <class ElemType, class option>
struct view_tester
{
	struct type
	{
		void operator() (mesh_t const& mesh)
		{
			std::cout << std::endl << "Elem type id: " << ElemType::id << std::endl;

			ElemType e(ElemType::coords_t::Random());
			std::cout << "coords: " << e.get_coords() << std::endl;

			auto cfv = constant_view(e);
			std::cout << cfv.get_dofs() << std::endl;
			
			auto ifv = isoparametric_view(e);
			std::cout << ifv.get_dofs() << std::endl;

			auto dcfv = dirac(cfv);
			std::cout << dcfv.get_dofs() << std::endl;
			
			auto difv = dirac(ifv);
			std::cout << difv.get_dofs() << std::endl;
		}
	};
	
};

int main(void)
{
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

	uMatrix elements(5, 1+4);
	elements <<
		quad_1_elem::id, 0, 1, 4, 3,
		quad_1_elem::id, 1, 2, 5, 4,
		quad_1_elem::id, 3, 4, 7, 6,
		tria_1_elem::id, 4, 5, 8, 0,
		tria_1_elem::id, 4, 8, 7, 0;

	mesh_t msh(nodes, elements);
	return 0;
}