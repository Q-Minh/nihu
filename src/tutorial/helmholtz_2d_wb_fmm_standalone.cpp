#include "core/field.hpp"
#include "core/function_space.hpp"
#include "fmm/helmholtz_exterior_solver.hpp"
#include "fmm/helmholtz_2d_field_point.hpp"
#include "fmm/helmholtz_2d_wb_fmm.hpp"
#include "interface/read_off_mesh.hpp"
#include "library/lib_element.hpp"
#include "util/type2tag.hpp"

#include <Eigen/IterativeLinearSolvers>
#include "fmm/GMRES.h"

// basic type parameter inputs
typedef NiHu::line_1_tag tag_t;
typedef double wave_number_t;

// computing the fmm type
typedef NiHu::fmm::helmholtz_2d_wb_fmm<wave_number_t> fmm_t;

// computing the fmbem type
typedef NiHu::tag2type<tag_t>::type elem_t;
typedef NiHu::field_view<elem_t, NiHu::field_option::constant> trial_field_t;
typedef NiHu::dirac_field<trial_field_t> test_field_t;
typedef elem_t::x_t location_t;

// computing the cluster tree type
typedef fmm_t::cluster_t cluster_t;
typedef NiHu::fmm::cluster_tree<cluster_t> cluster_tree_t;

typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;

void read_excitation(std::string fname, cvector_t &xct, double &k)
{
	std::ifstream ifs(fname);
	if (!ifs)
		throw std::runtime_error("Could not open inut file " + fname);
	ifs >> k;
	size_t N;
	ifs >> N;

	xct.resize(N, 1);
	for (size_t i = 0; i < N; ++i)
	{
		double r, im;
		ifs >> r >> im;
		xct(i, 0) = std::complex<double>(r, im);
	}
	ifs.close();
}

void export_response(std::string fname, cvector_t const &res, double k)
{
	std::ofstream ofs(fname);
	ofs << k << '\n';
	ofs << res.rows() << '\n';
	for (Eigen::Index i = 0; i < res.rows(); ++i)
		ofs << res(i, 0).real() << '\t' << res(i, 0).imag() << '\n';
	ofs.close();
}

int main(int argc, char *argv[])
{
	if (argc < 4)
	{
		std::cerr << "Use: " << argv[0] << " meshname fieldname xct_noext" << std::endl;
		return 1;
	}

	try
	{
		// read parameters
		std::string surf_mesh_name(argv[1]);
		std::string field_mesh_name(argv[2]);
		std::string xct_noext(argv[3]);
		std::string xct_name = xct_noext + ".xct";
		std::string surf_res_name = xct_noext + "_surf.res";
		std::string field_res_name = xct_noext + "_field.res";

		std::cout << "mesh: " << surf_mesh_name << std::endl;
		std::cout << "field: " << field_mesh_name << std::endl;
		std::cout << "xct name: " << xct_name << std::endl;
		std::cout << "surf res name: " << surf_res_name << std::endl;
		std::cout << "field res name: " << field_res_name << std::endl;

		// read mesh
		auto surf_mesh = NiHu::read_off_mesh(surf_mesh_name, tag_t());

		// read excitation
		double k;
		cvector_t q_surf, p_surf;
		read_excitation(xct_name, q_surf, k);

		std::cout << "wave number: " << k << std::endl;

		auto const &trial_space = NiHu::constant_view(surf_mesh);
		typedef std::decay<decltype(trial_space)>::type trial_space_t;

		NiHu::fmm::divide_num_nodes div(10);
		size_t far_field_quadrature_order = 6;

		// solve surface system
		{
			typedef NiHu::fmm::helmholtz_2d_wb_fmm<double> fmm_t;
			NiHu::fmm::helmholtz_exterior_solver<fmm_t, trial_space_t> solver(trial_space);
			solver.set_excitation(q_surf);
			solver.set_wave_number(k);
			p_surf = solver.solve(div, far_field_quadrature_order);
			export_response(surf_res_name, p_surf, k);
		}

		// read field
		auto field_mesh = NiHu::read_off_mesh(field_mesh_name, tag_t());

		// evaluate field pressure
		{
			auto const &test_space = NiHu::dirac(NiHu::constant_view(field_mesh));
			typedef std::decay<decltype(test_space)>::type test_space_t;

			NiHu::fmm::helmholtz_2d_field_point<test_space_t, trial_space_t> solver(test_space, trial_space);
			solver.set_qsurf(q_surf);
			solver.set_psurf(p_surf);
			solver.set_wave_number(k);
			solver.eval(div, far_field_quadrature_order);
			export_response(field_res_name, solver.get_response(), k);
		}
	}
	catch (std::exception const &e)
	{
		std::cerr << e.what() << std::endl;
	}
	catch (char const *str)
	{
		std::cerr << str << std::endl;
	}
	catch (...)
	{
		std::cerr << "Unhandled exception caught" << std::endl;
	}

	return 0;
}
