#include <iostream>
#include <fstream>
#include <sstream>

#include <omp.h>

#include "core/field.hpp"
#include "core/weighted_residual.hpp"
#include "interface/read_off_mesh.hpp"
#include "library/lib_element.hpp"
#include "library/quad_1_gauss_field.hpp"
#include "library/helmholtz.hpp"

#include<Eigen/IterativeLinearSolvers>
#include "fmm/GMRES.h"

#include <boost/math/constants/constants.hpp>
#include <boost/program_options.hpp>

#define GAUSS
#define BM
#define ITERATIVE

#ifdef BM
#define KHIE
#define HSIE
#endif

typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cVector;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> dVector;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;

template <class TestSpace, class TrialSpace>
cVector solve(TestSpace const &test_space, TrialSpace const &trial_space,
	double k, cMatrix const &qsurf, double tolerance, int restart, int &iterations)
{
	size_t nDof = trial_space.get_num_dofs();

#if defined(HSIE)
#if !defined(KHIE) 
	double alpha = 1.0;
#else
	std::complex<double> alpha(0., -1. / k);
#endif
#endif		

	// create integral operators
	auto Gop = NiHu::create_integral_operator(NiHu::helmholtz_3d_SLP_kernel<double>(k));
	auto Hop = NiHu::create_integral_operator(NiHu::helmholtz_3d_DLP_kernel<double>(k));
	auto Iop = NiHu::identity_integral_operator();

#ifdef HSIE
	auto Htop = NiHu::create_integral_operator(NiHu::helmholtz_3d_DLPt_kernel<double>(k));
	auto Dop = NiHu::create_integral_operator(NiHu::helmholtz_3d_HSP_kernel<double>(k));
#endif

	// compute rhs
	cVector rhs;

	{
		cMatrix Gs(nDof, nDof);
		Gs.setZero();

#ifdef KHIE			
		std::cout << "Integrating Gs" << std::endl;
		Gs << test_space * Gop[trial_space];
#endif

#ifdef HSIE
		std::cout << "Integrating Hts" << std::endl;
		Gs << test_space * (alpha * Htop)[trial_space];
		Gs << test_space * (alpha / 2. * Iop)[trial_space];
#endif			

		rhs = Gs * qsurf;
	}

	cVector psurf;

	{
		// create matrices
		cMatrix Hs(nDof, nDof);
		Hs.setZero();

#ifdef HSIE
		std::cout << "Integrating Ds" << std::endl;
		Hs << test_space * (alpha * Dop)[trial_space];
#endif

#ifdef KHIE
		std::cout << "Integrating Hs" << std::endl;
		Hs << test_space * (-.5 * Iop)[trial_space];
		Hs << test_space * Hop[trial_space];
#endif

		// solve linear system
		std::cout << "Solving linear system" << std::endl;
#ifdef ITERATIVE

		Eigen::GMRES<decltype(Hs), Eigen::IdentityPreconditioner > solver(Hs);
		solver.setTolerance(tolerance);
		solver.set_restart(restart);
		psurf = solver.solve(rhs);
		iterations = solver.iterations();

		std::cout << "k:               " << k << '\n';
		std::cout << "#iterations:     " << solver.iterations() << '\n';
		std::cout << "estimated error: " << solver.error() << std::endl;

#else
		psurf = Hs.colPivHouseholderQr().solve(rhs);
		std::cout << "k:               " << k << '\n';
		iterations = 1;
#endif
	}

	return psurf;
}


template <class TestSpace, class TrialSpace>
cVector postproc(TestSpace const &test_space, TrialSpace const &trial_space,
	double k, cMatrix const &psurf, cMatrix const &qsurf)
{
	size_t nDof = trial_space.get_num_dofs();
	size_t M = test_space.get_num_dofs();

	// create integral operators
	auto Gop = NiHu::create_integral_operator(NiHu::helmholtz_3d_SLP_kernel<double>(k));
	auto Hop = NiHu::create_integral_operator(NiHu::helmholtz_3d_DLP_kernel<double>(k));

	cMatrix Gf(M, nDof), Hf(M, nDof);
	Gf.setZero();
	Hf.setZero();

	std::cout << "Integrating Gf" << std::endl;
	Gf << test_space * Gop[trial_space];
	std::cout << "Integrating Hf" << std::endl;
	Hf << test_space * Hop[trial_space];

	return Hf * psurf - Gf * qsurf;
}


void read_excitation(std::string fname, cVector &xct, double &k)
{
	std::ifstream ifs(fname);
	if (!ifs)
		throw std::runtime_error("Could not open file for reading: " + fname);
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


void export_response(std::string fname, cVector const &res, double const &k, size_t iter = 1)
{
	std::ofstream ofs(fname);
	if (!ofs)
		throw std::runtime_error("Could not open file for writing: " + fname);
	ofs << k << '\n';
	ofs << res.rows() << '\n';
	for (Eigen::Index i = 0; i < res.rows(); ++i)
		ofs << res(i, 0).real() << '\t' << res(i, 0).imag() << '\n';
	ofs << iter << '\n';
	ofs.close();
}


int main(int argc, char **argv)
{
	using namespace boost::math::double_constants;
	using namespace boost::program_options;

	try
	{
		options_description desc("Options");
		desc.add_options()
			("help", "Help screen")
			("solve", "Solve surface system")
			("postprocess", "Compute field point pressure")
			("surface_mesh", value<std::string>(), "Surface mesh file name")
			("field_mesh", value<std::string>(), "Field point mesh file name")
			("surface_result", value<std::string>(), "Surface result file name")
			("field_result", value<std::string>(), "Field point result file name")
			("frequency", value<double>(), "Frequency [Hz]")
			("speed_of_sound", value<double>()->default_value(340), "Speed of sound [m/s]")
			("surface_velocity", value<double>()->default_value(1.0e-3), "Constant surface velocity excitation [m/s]")
			("density", value<double>()->default_value(1.3), "Density of air [kg/m3]")
			("tolerance", value<double>()->default_value(1e-8), "Tolerance of GMRES [-]")
			("restart", value<size_t>()->default_value(3000), "Restart parameter of GMRES [-]");

		variables_map vm;
		store(parse_command_line(argc, argv, desc), vm);
		notify(vm);

		if (vm.count("help"))
			std::cout << desc << '\n';

		// read parameters
		std::string surf_mesh_name = vm["surface_mesh"].as<std::string>();
		std::string surface_result_name = vm["surface_result"].as<std::string>();
		double freq = vm["frequency"].as<double>();
		double rho = vm["density"].as<double>();
		double c = vm["speed_of_sound"].as<double>();
		double v0 = vm["surface_velocity"].as<double>();
		size_t restart = vm["restart"].as<size_t>();
		double tolerance = vm["tolerance"].as<double>();

#ifdef GAUSS
		// read mesh file
		uMatrix elements;
		dMatrix nodes;
		NiHu::read_off_data(surf_mesh_name, nodes, elements, NiHu::quad_1_tag());

		// assemble field matrix
		size_t nElements = elements.rows();
		uMatrix fields(nElements, 1 + 4 + 4);
		for (size_t e = 0; e < nElements; ++e)
		{
			fields(e, 0) = NiHu::quad_1_gauss_field::id;
			for (size_t c = 0; c < 4; ++c)
				fields(e, c + 1) = elements(e, c + 1);
			for (size_t c = 0; c < 4; ++c)
				fields(e, c + 1 + 4) = 4 * e + c;
		}

		// create function space
		auto trial_space = NiHu::create_function_space(nodes, fields, NiHu::quad_1_gauss_field_tag());
#else
		auto mesh = NiHu::read_off_mesh(surf_mesh_name, NiHu::quad_1_tag(), NiHu::tria_1_tag());
		auto const &trial_space = NiHu::constant_view(mesh);
#endif
		auto const &test_space = NiHu::dirac(trial_space);

		// solve here
		double z0 = rho * c;
		double om = two_pi * freq;
		double k = om / c;

		std::complex<double> const J(0., 1.);

		cVector q_surf;
		q_surf.resize(trial_space.get_num_dofs());
		q_surf.setConstant(-J * k * z0 * v0);

		if (vm.count("solve") > 0)
		{
			std::cout << "mesh: " << surf_mesh_name << std::endl;
			std::cout << "result: " << surface_result_name << std::endl;
			std::cout << "freq: " << freq << std::endl;

			// solve
			int iters;
			cVector p_surf = solve(test_space, trial_space, k, q_surf, tolerance, restart, iters);

			// export ps
			export_response(surface_result_name, p_surf, k, iters);
		}

		if (vm.count("postprocess") > 0)
		{
			// read parameters
			std::string field_mesh_name = vm["field_mesh"].as<std::string>();
			std::string field_point_result_name = vm["field_result"].as<std::string>();

			std::cout << "mesh: " << surf_mesh_name << std::endl;
			std::cout << "field: " << field_mesh_name << std::endl;
			std::cout << "freq: " << freq << std::endl;

			// import ps
			cVector p_surf;
			double new_k;
			read_excitation(surface_result_name, p_surf, new_k);

			// field point pressure
			auto field = NiHu::read_off_mesh(field_mesh_name, NiHu::quad_1_tag());
			// create field function space
			auto const &field_sp = NiHu::dirac(NiHu::constant_view(field));

			// compute pf
			cVector p_field = postproc(field_sp, trial_space, k, q_surf, p_surf);

			// export pf
			export_response(field_point_result_name, p_field, k);
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
		std::cerr << "unhandled exception" << std::endl;
	}

	return 0;
}
