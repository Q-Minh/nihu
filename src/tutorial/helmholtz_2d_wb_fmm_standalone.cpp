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

// basic type parameter inputs

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
		typedef NiHu::line_1_tag tag_t;
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

		typedef double wave_number_t;
		typedef NiHu::fmm::helmholtz_2d_wb_fmm<wave_number_t> fmm_t;

		// solve surface system
		{
			NiHu::fmm::helmholtz_burton_miller_solver<fmm_t, trial_space_t> solver(trial_space);
			solver.set_excitation(q_surf);
			solver.set_wave_number(k);
			solver.set_far_field_order(far_field_quadrature_order);
			p_surf = solver.solve(div);
			export_response(surf_res_name, p_surf, k);
		}

		// read field
		auto field_mesh = NiHu::read_off_mesh(field_mesh_name, tag_t());

		// evaluate field pressure
		{
			auto const &test_space = NiHu::dirac(NiHu::constant_view(field_mesh));
			typedef std::decay<decltype(test_space)>::type test_space_t;

			auto cpu_t0 = NiHu::cpu_time::tic();
			auto wc_t0 = NiHu::wc_time::tic();

			NiHu::fmm::helmholtz_field_point<fmm_t, test_space_t, trial_space_t> evaluator(test_space, trial_space);
			evaluator.set_qsurf(q_surf);
			evaluator.set_psurf(p_surf);
			evaluator.set_wave_number(k);
			evaluator.set_far_field_order(far_field_quadrature_order);
			evaluator.eval(div);
			export_response(field_res_name, evaluator.get_response(), k);

			double cpu_t = NiHu::cpu_time::toc(cpu_t0);
			double wc_t = NiHu::wc_time::toc(wc_t0);

			std::cout << "Elapsed cpu time: " << cpu_t << " s" << std::endl;
			std::cout << "Elapsed wall clock time: " << wc_t << " s"<< std::endl;
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
