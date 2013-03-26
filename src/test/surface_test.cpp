#include "..\bem\double_integral.hpp"

#include "..\bem\field.hpp"

typedef quad_1_elem elem_t;
typedef field<elem_t, isoparametric_field> field_t;
typedef green_G_kernel kernel_t;

int main(void)
{
	elem_t::coords_t test_coords;
	test_coords <<
		0.0, 0.0, 0.0,
		1.0, 0.0, 0.0,
		1.0, 1.0, 0.0,
		0.0, 1.0, 0.0;

	elem_t::coords_t trial_coords;
	trial_coords <<
		0.0, 0.0, 1.0,
		1.0, 0.0, 1.0,
		1.0, 1.0, 1.0,
		0.0, 1.0, 1.0;

	elem_t test_elem(test_coords), trial_elem(trial_coords);
	field_t test_field(test_elem), trial_field(trial_elem);

	kernel_t::set_wave_number(1.0);

	double_integral<kernel_t, field_t, field_t> di;

	std::cout << di.eval(test_field, trial_field);

	return 0;
}
