#include <iostream>
#include "core/field.hpp"
#include "tmp/sequence.hpp"
#include "tmp/control.hpp"
#include "tmp/vector.hpp"

template <class ElemType, class option>
struct view_tester
{
	struct type
	{
		void operator() (void)
		{
			std::cout << std::endl << "Elem type id: " << ElemType::id << std::endl;

			ElemType e(ElemType::coords_t::Random());

			auto fv = create_field_view(e, option());
			std::cout << fv.get_dofs() << std::endl;
			
			auto dfv = dirac(fv);
			std::cout << dfv.get_dofs() << std::endl;
		}
	};
};

int main(void)
{
	typedef tmp::vector<tria_1_elem, quad_1_elem, tria_2_elem, quad_2_elem> elem_vector;
	typedef tmp::vector<field_option::constant, field_option::isoparametric> option_vector;

	tmp::d_call_each<
		elem_vector,
		option_vector,
		view_tester<tmp::_1, tmp::_2>
	>();

	return 0;
}

