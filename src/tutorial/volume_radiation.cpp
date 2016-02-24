#include "core/weighted_residual.hpp"
#include "interface/read_off_mesh.hpp"
#include "library/helmholtz_kernel.hpp"

#include <sstream>
#include <string>

typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cMatrix;

typedef volume_element<quad_1_shape_set, double> vol_quad_elem;
struct vol_quad_tag {};
template <> struct tag2element<vol_quad_tag> { typedef vol_quad_elem type; };


typedef helmholtz_2d_double_kernel<std::complex<double>, 1, 1> Q11Kernel;
typedef helmholtz_2d_double_kernel<std::complex<double>, 1, 2> Q12Kernel;
typedef helmholtz_2d_double_kernel<std::complex<double>, 2, 2> Q22Kernel;

// main(volume_mesh.off, field_mesh.off, volume_excitation.ascii, wave_number_real, wave_number_imag
int main(int argc, char *argv[])
{
	// import arguments
	if (argc != 6)
	{
		std::stringstream msg;
		msg << "Usage: " << argv[0] << " volume_mesh.off field_mesh.off volume_excitation.ascii wave_number_real wave_number_imag";
		throw std::runtime_error(msg.str());
	}
	char const *volume_mesh_file = argv[1];
	char const *field_mesh_file = argv[2];
	char const *excitation_mesh_file = argv[3];
	std::complex<double> k(std::stod(argv[4]), std::stod(argv[5]));
	
	// read meshes
	auto domain_mesh = read_off_mesh("volume_mesh.off", vol_quad_tag());
	auto field_mesh = read_off_mesh("surface_mesh.off", line_1_tag());

	// create function spaces
	auto const &trial_space = create_function_space_view(domain_mesh, field_option::constant());
	auto const &test_space = dirac(constant_view(field_mesh));
	auto nField = trial_space.get_num_dofs();
	auto nVol = test_space.get_num_dofs();

	// read excitation
	cMatrix xct(nVol, 2);
	std::ifstream is(excitation_mesh_file);
	for (int i = 0; i < nVol; ++i)
	{
		double ur, ui, vr, vi;
		is >> ur >> ui >> vr >> vi;
		xct.row(i) << std::complex<double>(ur, ui), std::complex<double>(vr, vi);
	}
	is.close();

	auto G11_op = create_integral_operator(Q11Kernel(k));
	auto G12_op = create_integral_operator(Q12Kernel(k));
	auto G22_op = create_integral_operator(Q22Kernel(k));

	cMatrix G(nField, nVol);
	cMatrix resp(nField, 1);
	resp.setZero();

	G.setZero();
	G << test_space * G11_op[trial_space];
	resp += G * xct.col(0).array().square().matrix();
	G.setZero();
	G << test_space * G12_op[trial_space];
	resp += 2. * G * (xct.col(0).array() * xct.col(1).array()).matrix();
	G.setZero();
	G << test_space * G22_op[trial_space];
	resp += G * xct.col(1).array().square().matrix();

	std::cout << resp;

	return 0;
}

