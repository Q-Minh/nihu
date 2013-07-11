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

			typedef field_view<ElemType, option> fv_t;
			std::cout << static_cast<fv_t const &>(e).get_dofs() << std::endl;
			std::cout << field_traits<fv_t>::is_dirac << std::endl;

			typedef dirac_field<fv_t> dv_t;
			std::cout << static_cast<dv_t const &>(e).get_dofs() << std::endl;
			std::cout << field_traits<dv_t>::is_dirac << std::endl;
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

