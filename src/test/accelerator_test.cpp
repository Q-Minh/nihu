
#include "../bem/element.hpp"
#include "../bem/kernel.hpp"
#include "../bem/accelerator.hpp"

typedef accelerator<quad_1_elem, green_G_kernel::input_t> accelerator_t;

int main(void)
{
	accelerator_t a;
	return 0;
}
