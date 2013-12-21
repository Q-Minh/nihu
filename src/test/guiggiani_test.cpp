#include <library/guiggiani_1992.hpp>
#include <core/element.hpp>
#include <core/field.hpp>
#include <library/laplace_kernel.hpp>

int main(void)
{
	typedef quad_1_elem elem_t;
	typedef field_view<elem_t, field_option::constant> field_t;
	typedef laplace_3d_HSP_kernel kernel_t;
	typedef elem_t::xi_t xi_t;

	typedef guiggiani_hypersingular_collocation<kernel_t, field_t> guiggiani_t;

	elem_t elem(elem_t::coords_t::Random());
	xi_t xi(xi_t::Constant(.3));
	double theta(1.0);

	std::cout << "N0: " << guiggiani_t::N0(xi) << std::endl;
	std::cout << "N1: " << guiggiani_t::N1(xi, theta) << std::endl;
	std::cout << "A vector: " << guiggiani_t::A_vector(elem, xi, theta).transpose() << std::endl;
	std::cout << "B vector: " << guiggiani_t::B_vector(elem, xi, theta).transpose() << std::endl;

	return 0;
}
