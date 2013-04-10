#include <iostream>
#include "../bem/singularity_check.hpp"

typedef quad_1_elem elem_t;
typedef field<elem_t, constant_field, dirac_field> test_field_t;
typedef field<elem_t, constant_field, function_field> trial_field_t;

class kernel {};

int main(void)
{
	elem_t test_elem = elem_t(elem_t::coords_t());
	test_field_t test_field(test_elem);
	trial_field_t trial_field(test_elem);

	std::cout << (singularity_check<kernel, test_field_t, trial_field_t>::eval(test_field, trial_field) == FACE_MATCH) << std::endl;

	return 0;
}
