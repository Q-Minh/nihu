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

			auto fv = create_field_view(e, option());
			std::cout << fv.get_dofs() << std::endl;
//			std::cout << "fv_t::is_dirac = " << std::boolalpha << field_traits<decltype(cfv)>::is_dirac << std::endl;

			auto dv = dirac(fv);
			std::cout << dv.get_dofs() << std::endl;
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

	tmp::call_each<
		elem_vector,
		view_tester<tmp::_1, field_option::isoparametric>
	>();

	return 0;
}

