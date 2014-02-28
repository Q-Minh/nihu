#include "library/elastostatics_kernel.hpp"
#include "library/lib_element.hpp"

int main(void)
{
    elastostatics_3d_U_kernel U(.33);

    quad_1_elem::coords_t coords1, coords2;
    coords1 <<
		0.0, 1.0, 1.0, 0.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;
    coords2 <<
		2.0, 3.0, 3.0, 1.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;

    quad_1_elem elem1(coords1), elem2(coords2);

    elastostatics_3d_U_kernel::test_input_t test_input(elem1, quad_1_elem::xi_t::Zero());
    elastostatics_3d_U_kernel::trial_input_t trial_input(elem2, quad_1_elem::xi_t::Zero());

    std::cout << U(test_input, trial_input);
}
