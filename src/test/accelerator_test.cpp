
#include "../bem/element.hpp"
#include "../bem/kernel.hpp"
#include "../bem/accelerator.hpp"

#include "../tmp/vector.hpp"
#include "../tmp/control.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_type_vector_t;
typedef accelerator<elem_type_vector_t, green_G_kernel::input_t> accelerator_t;

int main(void)
{
	accelerator_t a;
	for (unsigned i = 0; i < 2000; ++i)
	{
		quad_1_elem elem(quad_1_elem::coords_t::Random());
		a.add_elem(elem);

		tria_1_elem telem(tria_1_elem::coords_t::Random());
		a.add_elem(telem);
	}

	std::for_each(a.elem_begin<quad_1_elem>(), a.elem_end<quad_1_elem>(),
		[] (accelerator_t::elem_pool_t<quad_1_elem>::type const &ep) {
			std::for_each(ep[0].begin(), ep[0].end(),
				[] (green_G_kernel::input_t const & i) {
					std::cout << i.get_jacobian() << std::endl;
			});
	});

	return 0;
}
