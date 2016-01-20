// #include "core/element.hpp"
#include "core/weighted_residual.hpp"
#include "interface/read_off_mesh.hpp"
// #include "library/helmholtz_kernel.hpp"


typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cMatrix;


typedef volume_element<quad_1_shape_set, double> vol_quad_elem;

struct vol_quad_tag {};

template <>
struct tag2element<vol_quad_tag>
{
	typedef vol_quad_elem type;
};



struct MyKernel
{
	typedef Eigen::Matrix<std::complex<double>, 1, 3> return_type;
	return_type operator()(location_input_2d const &x, location_input_2d const &y, empty_data const &)
	{
		return return_type::Constant(1./(x.get_x()-y.get_x()).norm());
	}
};


class MyVeryKernel;

template <>
struct kernel_traits<MyVeryKernel>
{
	typedef location_input_2d test_input_t;
	typedef location_input_2d trial_input_t;
	typedef collect<empty_data> data_t;
	typedef single_brick_wall<MyKernel>::type output_t;
	typedef asymptotic::inverse<1> far_field_behaviour_t;
	typedef gauss_family_tag quadrature_family_t;
	static bool const is_symmetric = true;
	static bool const is_singular = true;
	enum { result_rows = 1, result_cols = 3 };
};

class MyVeryKernel : public kernel_base<MyVeryKernel>
{
};


int main(void)
{
	auto domain_mesh = read_off_mesh("volume_mesh.off", vol_quad_tag());
	auto surface_mesh = read_off_mesh("surface_mesh.off", line_1_tag());

	auto elem = *domain_mesh.begin<vol_quad_elem>();

	auto const &trial_space = create_function_space_view(domain_mesh, field_option::constant(), _3d());
	auto const &test_space = dirac(constant_view(surface_mesh));

	double k = 1.0;

	auto G_op = create_integral_operator(MyVeryKernel());

	cMatrix G(test_space.get_num_dofs(), trial_space.get_num_dofs());
	G.setZero();
	G << test_space * G_op[trial_space];

	std::cout << G;

	return 0;
}

