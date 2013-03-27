#include "..\bem\double_integral.hpp"

#include "..\bem\field.hpp"

typedef tria_1_elem test_elem_t;
typedef field<test_elem_t, constant_field> test_field_t;

typedef quad_1_elem trial_elem_t;
typedef field<trial_elem_t, constant_field> trial_field_t;

typedef green_G_kernel kernel_t;

int main(void)
{
	test_elem_t::coords_t test_coords;
	test_coords <<
		0.0, 0.0, 0.0,
		1.0, 0.0, 0.0,
		1.0, 1.0, 0.0;

	trial_elem_t::coords_t trial_coords;
	trial_coords <<
		0.0, 0.0, 1.0,
		1.0, 0.0, 1.0,
		1.0, 1.0, 1.0,
		0.0, 1.0, 1.0;

	test_elem_t test_elem(test_coords);
	trial_elem_t trial_elem(trial_coords);
	test_field_t test_field(test_elem);
	trial_field_t trial_field(trial_elem);

	kernel_t::set_wave_number(1.0);

	double_integral<kernel_t, test_field_t, trial_field_t> di;

	std::cout << di.eval(test_field, trial_field);

	return 0;
}
