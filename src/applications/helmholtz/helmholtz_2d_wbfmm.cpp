#include "core/field.hpp"
#include "core/function_space.hpp"
#include "fmm/helmholtz_burton_miller_solver.hpp"
#include "fmm/helmholtz_field_point.hpp"
#include "fmm/helmholtz_2d_wb_fmm.hpp"
#include "interface/read_off_mesh.hpp"
#include "library/lib_element.hpp"
#include "util/type2tag.hpp"
#include "util/timer.h"

#include <Eigen/IterativeLinearSolvers>
#include "fmm/GMRES.h"

#include <boost/program_options.hpp>


// basic type parameter inputs

typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;
typedef double wave_number_t;

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
	using namespace boost::program_options;

	try
	{
		options_description desc("Options");
		desc.add_options()
			("help", "Help screen")
			("solve", "Solve surface system")
			("postprocess", "Compute field point pressure")
			("surface_mesh", value<std::string>(), "Surface mesh file name")
			("surface_excitation", value<std::string>(), "Surface normal derivative exc.file name")
			("field_mesh", value<std::string>(), "Field point mesh file name")
			("surface_result", value<std::string>(), "Surface result file name")
			("field_result", value<std::string>(), "Field point result file name")
			("tolerance", value<double>()->default_value(1e-8), "Tolerance of GMRES [-]")
			("restart", value<size_t>()->default_value(3000), "Restart parameter of GMRES [-]")
			("num_leaf_nodes", value<size_t>()->default_value(10), "Number of leaf level nodes [-]")
			("far_field_order", value<size_t>()->default_value(6), "Order of far field quadrature [-]");

		variables_map vm;
		store(parse_command_line(argc, argv, desc), vm);
		notify(vm);

		if (vm.count("help")) {
			std::cout << desc << '\n';
			return 0;
		}

		// read parameters
		std::string surface_result_name = vm["surface_result"].as<std::string>();
		size_t restart = vm["restart"].as<size_t>();
		double tolerance = vm["tolerance"].as<double>();
		size_t far_field_quadrature_order = vm["far_field_order"].as<size_t>();
		size_t num_leaf_nodes = vm["num_leaf_nodes"].as<size_t>();

		// read parameters
		std::string surf_mesh_name = vm["surface_mesh"].as<std::string>();
		std::string surf_xct_name = vm["surface_excitation"].as<std::string>();
		std::string surf_res_name = vm["surface_result"].as<std::string>();

		// read mesh
		typedef NiHu::line_1_tag tag_t;
		auto surf_mesh = NiHu::read_off_mesh(surf_mesh_name, tag_t());

		auto const &trial_space = NiHu::constant_view(surf_mesh);
		typedef std::decay<decltype(trial_space)>::type trial_space_t;

		NiHu::fmm::divide_num_nodes div(num_leaf_nodes);

		typedef NiHu::fmm::helmholtz_2d_wb_fmm<wave_number_t> fmm_t;

		cvector_t q_surf, p_surf;
		wave_number_t xct_k;
		read_excitation(surf_xct_name, q_surf, xct_k);

		if (vm.count("solve") > 0)
		{
			NiHu::fmm::helmholtz_burton_miller_solver<fmm_t, trial_space_t> solver(trial_space);
			solver.set_excitation(q_surf);
			solver.set_wave_number(xct_k);
			solver.set_far_field_order(far_field_quadrature_order);
			solver.set_restart(restart);
			solver.set_tolerance(tolerance);
			p_surf = solver.solve(div);
			export_response(surf_res_name, p_surf, xct_k);
		}

		if (vm.count("postprocess") > 0)
		{
			std::string field_mesh_name = vm["field_mesh"].as<std::string>();
			std::string field_point_result_name = vm["field_result"].as<std::string>();

			// read surface result
			double surf_k;
			read_excitation(surface_result_name, p_surf, surf_k);

			if (std::abs(surf_k / xct_k - 1) > 1e-8)
				throw "Excitation and surface response wave number mismatch !!!";

			// read field
			auto field_mesh = NiHu::read_off_mesh(field_mesh_name, tag_t());

			auto const &test_space = NiHu::dirac(NiHu::constant_view(field_mesh));
			typedef std::decay<decltype(test_space)>::type test_space_t;

			auto cpu_t0 = NiHu::cpu_time::tic();
			auto wc_t0 = NiHu::wc_time::tic();

			NiHu::fmm::helmholtz_field_point<fmm_t, test_space_t, trial_space_t> field_bir(test_space, trial_space);
			field_bir.set_qsurf(q_surf);
			field_bir.set_psurf(p_surf);
			field_bir.set_wave_number(surf_k);
			field_bir.set_far_field_order(far_field_quadrature_order);
			field_bir.eval(div);
			export_response(field_point_result_name, field_bir.get_response(), surf_k);
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
