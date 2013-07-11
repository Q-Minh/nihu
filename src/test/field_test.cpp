#include <iostream>
#include "../bem/field.hpp"
#include "../tmp/sequence.hpp"
#include "../tmp/control.hpp"

template <class ElemType, class option>
struct view_tester
{
	struct type
	{
		void operator() (void)
		{
			std::cout << std::endl << "Elem type id: " << ElemType::id << std::endl;

			ElemType e(ElemType::coords_t::Random());
			std::cout << "coords: " << e.get_coords() << std::endl;

			auto cfv = constant_view(e);
			std::cout << cfv.get_dofs() << std::endl;
//			std::cout << "fv_t::is_dirac = " << std::boolalpha << field_traits<decltype(cfv)>::is_dirac << std::endl;

			auto dfv = dirac(cfv);
			std::cout << dfv.get_dofs() << std::endl;
//			std::cout << "dv_t::is_dirac = " << std::boolalpha << field_traits<dv_t>::is_dirac << std::endl;
		}
	};
};

int main(void)
{
	typedef tmp::vector<tria_1_elem, quad_1_elem, tria_2_elem, quad_2_elem> elem_vector;

	tmp::call_each<
		elem_vector,
		view_tester<tmp::_1, field_option::constant>
	>();

	return 0;
}

