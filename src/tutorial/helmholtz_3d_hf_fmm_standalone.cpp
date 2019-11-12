#include "core/field.hpp"
#include "core/function_space.hpp"
#include "fmm/divide.hpp"
#include "fmm/helmholtz_3d_hf_fmm.hpp"
#include "fmm/helmholtz_burton_miller_solver.hpp"
#include "fmm/helmholtz_field_point.hpp"
#include "interface/read_off_mesh.hpp"
#include "library/lib_element.hpp"
#include "library/quad_1_gauss_field.hpp"

#include <boost/math/constants/constants.hpp>
#include <boost/program_options.hpp>

#include <cstdlib>

//#define GAUSS

typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cVector;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> dVector;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;

#ifdef GAUSS
void read_off_data(std::string const &fname, dMatrix &nodes, uMatrix &elements)
{
	// open mesh file for reading
	std::ifstream is(fname);
	if (!is)
		throw std::runtime_error("Error reading mesh file");

	// read header from file (first row is 'OFF')
	std::string header;
	if (!(is >> header) || header != "OFF")
		throw std::runtime_error("Possibly invalid off file");

	// read number of nodes and number of elements, nEdges is dropped
	unsigned nNodes, nElements, nEdges;
	if (!(is >> nNodes >> nElements >> nEdges))
		throw std::runtime_error("Error reading number of mesh entries");

	// read nodes
	nodes.resize(nNodes, 3);
	for (unsigned i = 0; i < nNodes; ++i)
		if (!(is >> nodes(i, 0) >> nodes(i, 1) >> nodes(i, 2)))
			throw std::runtime_error("Error reading mesh nodes");

	// read elements
	elements.resize(nElements, 5);
	for (unsigned i = 0; i < nElements; ++i)
	{
		unsigned nvert;
		if (!(is >> nvert))
			throw std::runtime_error("Error reading mesh elements");
		for (unsigned c = 0; c < nvert; ++c)
			if (!(is >> elements(i, c + 1)))
				throw std::runtime_error("Error reading mesh elements");
		elements(i, 0) = NiHu::quad_1_elem::id;
	}

	is.close();
}
#endif


// basic type parameter inputs
typedef double wave_number_t;

// computing the fmbem type
#ifdef GAUSS
typedef NiHu::quad_1_gauss_field trial_field_t;
#else
typedef NiHu::field_view<NiHu::quad_1_elem, NiHu::field_option::constant> trial_field_t;
#endif
typedef trial_field_t::elem_t elem_t;
typedef elem_t::x_t location_t;


typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;

void read_excitation(std::string fname, cvector_t &xct, wave_number_t &k)
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

void export_response(std::string fname, cvector_t const &res, wave_number_t const &k, size_t iter = 1)
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


int main(int argc, char *argv[])
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
			("density", value<double>()->default_value(1.3), "Density of air [kg/m3]");

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

#ifdef GAUSS
		// read mesh file
		uMatrix elements;
		dMatrix nodes;
		read_off_data(surf_mesh_name, nodes, elements);

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
		auto mesh = NiHu::read_off_mesh(surf_mesh_name, NiHu::quad_1_tag());
		auto const &trial_space = NiHu::constant_view(mesh);
		typedef std::decay<decltype(trial_space)>::type trial_space_t;
#endif

		std::complex<double> const J(0., 1.);

		// generate excitation
		double z0 = rho * c;
		double om = two_pi * freq;
		wave_number_t k = om / c;

		cvector_t q_surf;
		q_surf.resize(trial_space.get_num_dofs());
		q_surf.setConstant(-J * k * z0 * v0);

		typedef NiHu::fmm::helmholtz_3d_hf_fmm<wave_number_t> fmm_t;

		double leaf_diameter = 1. / std::real(k);
		size_t far_field_quadrature_order = 6;

		if (vm.count("solve") > 0)
		{
			std::cout << "mesh: " << surf_mesh_name << std::endl;
			std::cout << "result: " << surface_result_name << std::endl;
			std::cout << "freq: " << freq << std::endl;

			auto solver = NiHu::fmm::create_helmholtz_burton_miller_solver(NiHu::type2tag<fmm_t>(), trial_space);
			solver.set_wave_number(k);
			solver.set_excitation(q_surf);
			cvector_t p_surf = solver.solve(NiHu::fmm::divide_diameter(leaf_diameter), far_field_quadrature_order);

			// export ps
			export_response(surface_result_name, p_surf, k, solver.get_iterations());
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
			cvector_t p_surf;
			wave_number_t new_k;
			read_excitation(surface_result_name, p_surf, new_k);

			// field point pressure
			auto field = NiHu::read_off_mesh(field_mesh_name, NiHu::quad_1_tag());
			auto const &test_space = NiHu::dirac(NiHu::constant_view(field));
			typedef std::decay<decltype(test_space)>::type test_space_t;

			auto field_bie = NiHu::fmm::create_helmholtz_field_point(
				NiHu::type2tag<fmm_t>::type(), test_space, trial_space);
			field_bie.set_wave_number(k);
			field_bie.set_psurf(p_surf);
			field_bie.set_qsurf(q_surf);
			auto p_field = field_bie.eval(NiHu::fmm::divide_diameter(leaf_diameter), far_field_quadrature_order);

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
