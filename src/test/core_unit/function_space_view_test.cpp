#include "tmp/vector.hpp"
#include "bem/function_space.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_type_vector_t;
typedef mesh<elem_type_vector_t> mesh_t;

template <class ElemType, class Option>
struct view_tester
{
	struct type
	{
		void operator() (mesh_t const& mesh)
		{
			auto& fsp_view  = create_function_space_view(mesh, Option());
			std::cout << "# DOFS: " << fsp_view.get_num_dofs() << std::endl;
			
			typedef field_view<ElemType, Option> field_view_t;
			
			auto b = fsp_view.template field_begin<field_view_t>();
			auto e = fsp_view.template field_end<field_view_t>();
			
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

typedef tmp::vector<field_option::constant, field_option::isoparametric> option_vector;

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
	
	tmp::d_call_each<
		elem_type_vector_t,
		option_vector,
		view_tester<tmp::_1, tmp::_2>,
		mesh_t const &
	>(msh);
	
	return 0;
}