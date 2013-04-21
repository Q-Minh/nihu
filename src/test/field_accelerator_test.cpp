#include "../bem/field_type_accelerator.hpp"
#include "../bem/gaussian_quadrature.hpp"

typedef quad_1_elem elem_t;
typedef field<elem_t, isoparametric_field, function_field> field_t;

#include <iostream>

int main(void)
{
	typedef gauss_family_tag gauss;
	
	field_type_accelerator<field_t, gauss> fta(5);
	for (auto it = fta.cbegin(); it != fta.cend(); ++it)
	{
		std::cout << it->get_quadrature_elem().get_xi().transpose() << std::endl;
		std::cout << it->get_shape().transpose() << std::endl;
	}
	return 0;
}
