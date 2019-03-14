#include "cluster_tree.hpp"
#include "divide.hpp"
#include "elem_center_iterator.hpp"
#include "fmm_matrix.hpp"
#include "helmholtz_3d_hf_fmm.hpp"
#include "matrix_free.hpp"
#include "p2p_integral.hpp"
#include "p2x_indexed.hpp"
#include "p2x_integral.hpp"
#include "x2p_indexed.hpp"
#include "x2p_integral.hpp"

#include "core/field.hpp"
#include "core/function_space.hpp"
#include "interface/read_off_mesh.hpp"
#include "library/lib_element.hpp"
#include "library/quad_1_gauss_field.hpp"

#include <Eigen/IterativeLinearSolvers>
#include "GMRES.h"

#ifndef NUM_PROCESSORS
#define NUM_PROCESSORS 1
#endif

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

// computing the fmm type
typedef fmm::helmholtz_3d_hf_fmm<wave_number_t> fmm_t;

// computing the fmbem type
#ifdef GAUSS
typedef NiHu::quad_1_gauss_field trial_field_t;
#else
typedef NiHu::field_view<NiHu::quad_1_elem, NiHu::field_option::constant> trial_field_t;
#endif
typedef trial_field_t::elem_t elem_t;
typedef elem_t::x_t location_t;

// computing the cluster tree type
typedef fmm_t::cluster_t cluster_t;
typedef fmm::cluster_tree<cluster_t> cluster_tree_t;

typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;

void read_excitation(std::string fname, cvector_t &xct, double &k)
{
	std::ifstream ifs(fname);
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

void export_response(std::string fname, cvector_t const &res, double k, int iter = 1)
{
	std::ofstream ofs(fname);
	ofs << k << '\n';
	ofs << res.rows() << '\n';
	for (size_t i = 0; i < res.rows(); ++i)
		ofs << res(i, 0).real() << '\t' << res(i, 0).imag() << '\n';
	ofs << iter << '\n';
	ofs.close();
}


int main(int argc, char *argv[])
{
	try
	{


		if (argc < 4)
		{
			std::cerr << "Use: " << argv[0] << " meshname fieldname pattern fstart" << std::endl;
			return 1;
		}

		// read parameters
		std::string surf_mesh_name(argv[1]);
		std::string field_mesh_name(argv[2]);
		std::string pattern(argv[3]);
		double freq = std::atof(argv[4]);

		std::cout << "mesh: " << surf_mesh_name << std::endl;
		std::cout << "field: " << field_mesh_name << std::endl;
		std::cout << "pattern: " << pattern << std::endl;

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
#endif


		// generate excitation
		double rho = 1.3;
		double c = 340.0;
		double v0 = 1e-3;
		double z0 = rho * c;

		std::complex<double> const J(0., 1.);

		// loop over frequencies
		double om = 2. * M_PI * freq;
		double k = om / c;

		cvector_t xct, p_surf;
		xct.resize(trial_space.get_num_dofs());
		xct.setConstant(-std::complex<double>(0, 1.0) * k * z0 * v0);
		p_surf.resize(trial_space.get_num_dofs());
		p_surf.setConstant(1.0);

		std::cout << "wave number: " << k << std::endl;
		// create kernel, fmm and fmbem objects
		fmm_t fmm(k);

		// get x2x fmm operators
		auto m2m_op = fmm.create_m2m();
		auto l2l_op = fmm.create_l2l();
		auto m2l_op = fmm.create_m2l();

		// p2x operators and derivatives
		auto p2m_op_0 = fmm.create_p2m<0>();
		auto p2l_op_0 = fmm.create_p2l<0>();
		auto p2m_op_1 = fmm.create_p2m<1>();
		auto p2l_op_1 = fmm.create_p2l<1>();

		// x2p operators and derivatives
		auto m2p_op_0 = fmm.create_m2p<0>();
		auto l2p_op_0 = fmm.create_l2p<0>();
		auto m2p_op_1 = fmm.create_m2p<1>();
		auto l2p_op_1 = fmm.create_l2p<1>();

		// p2p operators and derivatives
		auto p2p_op_00 = fmm.create_p2p<0, 0>();
		auto p2p_op_01 = fmm.create_p2p<0, 1>();
		auto p2p_op_10 = fmm.create_p2p<1, 0>();
		auto p2p_op_11 = fmm.create_p2p<1, 1>();

		// integrate operators over fields
		size_t quadrature_order = 10;


	// evaluate field pressure
	{
		typedef NiHu::dirac_field<NiHu::field_view<NiHu::quad_1_elem, NiHu::field_option::constant> > test_field_t;
		auto field_mesh = NiHu::read_off_mesh(field_mesh_name, NiHu::quad_1_tag());
		auto const &test_space = NiHu::dirac(NiHu::constant_view(field_mesh));

		double D = 2.8;
		size_t depth = unsigned(1 + std::log2(k*D));
#ifdef GAUSS
		cluster_tree_t tree(
			fmm::create_field_center_iterator(trial_space.field_begin<trial_field_t>()),
			fmm::create_field_center_iterator(trial_space.field_end<trial_field_t>()),
			fmm::create_elem_center_iterator(field_mesh.begin<elem_t>()),
			fmm::create_elem_center_iterator(field_mesh.end<elem_t>()),
			fmm::divide_depth(depth));
#else
		cluster_tree_t tree(
			fmm::create_elem_center_iterator(mesh.begin<elem_t>()),
			fmm::create_elem_center_iterator(mesh.end<elem_t>()),
			fmm::create_elem_center_iterator(field_mesh.begin<elem_t>()),
			fmm::create_elem_center_iterator(field_mesh.end<elem_t>()),
			fmm::divide_depth(depth));
#endif

		std::cout << tree << std::endl;


#if 0
		// solve surface system
		{
			typedef NiHu::dirac_field<trial_field_t> test_field_t;

			auto const &test_space = NiHu::dirac(trial_space);

			double D = 2.5;
			size_t depth = unsigned(std::log2(k*D));

#ifdef GAUSS
			cluster_tree_t tree(
				fmm::create_field_center_iterator(trial_space.field_begin<trial_field_t>()),
				fmm::create_field_center_iterator(trial_space.field_end<trial_field_t>()),
				fmm::create_field_center_iterator(test_space.field_begin<test_field_t>()),
				fmm::create_field_center_iterator(test_space.field_end<test_field_t>()),
				fmm::divide_depth(depth));
#else
			cluster_tree_t tree(
				fmm::create_elem_center_iterator(mesh.begin<elem_t>()),
				fmm::create_elem_center_iterator(mesh.end<elem_t>()),
				fmm::create_elem_center_iterator(mesh.begin<elem_t>()),
				fmm::create_elem_center_iterator(mesh.end<elem_t>()),
				fmm::divide_depth(depth));
#endif

			std::cout << tree << std::endl;

			// create interaction lists
			fmm::interaction_lists lists(tree);

			// x2p operators
			fmm::x2p_integral<decltype(m2p_op_0), test_field_t> im2p_op_0(m2p_op_0, quadrature_order);
			fmm::x2p_integral<decltype(l2p_op_0), test_field_t> il2p_op_0(l2p_op_0, quadrature_order);
			fmm::x2p_integral<decltype(m2p_op_1), test_field_t> im2p_op_1(m2p_op_1, quadrature_order);
			fmm::x2p_integral<decltype(l2p_op_1), test_field_t> il2p_op_1(l2p_op_1, quadrature_order);

			// p2p operators
			fmm::p2p_integral<decltype(p2p_op_00), test_field_t, trial_field_t> ip2p_op_00(p2p_op_00, true);
			fmm::p2p_integral<decltype(p2p_op_01), test_field_t, trial_field_t> ip2p_op_01(p2p_op_01, true);
			fmm::p2p_integral<decltype(p2p_op_10), test_field_t, trial_field_t> ip2p_op_10(p2p_op_10, true);
			fmm::p2p_integral<decltype(p2p_op_11), test_field_t, trial_field_t> ip2p_op_11(p2p_op_11, true);

			// create indexed fmbem operators
			std::complex<double> alpha(0.0, -1.0 / k);

			auto m2p_bm = fmm::create_x2p_indexed(im2p_op_0 + alpha * im2p_op_1, test_space.field_begin<test_field_t>(), test_space.field_end<test_field_t>());
			auto l2p_bm = fmm::create_x2p_indexed(il2p_op_0 + alpha * il2p_op_1, test_space.field_begin<test_field_t>(), test_space.field_end<test_field_t>());

			auto p2p_0 = fmm::create_p2x_indexed(
				fmm::create_x2p_indexed(ip2p_op_00 + alpha * ip2p_op_10,
					test_space.field_begin<test_field_t>(),
					test_space.field_end<test_field_t>()),
				trial_space.field_begin<trial_field_t>(),
				trial_space.field_end<trial_field_t>()
			);
			auto p2p_1 = fmm::create_p2x_indexed(
				fmm::create_x2p_indexed(ip2p_op_01 + alpha * ip2p_op_11,
					test_space.field_begin<test_field_t>(),
					test_space.field_end<test_field_t>()),
				trial_space.field_begin<trial_field_t>(),
				trial_space.field_end<trial_field_t>()
			);

			fmm.init_level_data(tree, 3.0);
			for (size_t c = 0; c < tree.get_n_clusters(); ++c)
				tree[c].set_p_level_data(&fmm.get_level_data(tree[c].get_level()));

			std::cout << "Starting assembling matrices " << std::endl;

			// create matrix objects
			auto slp_matrix = fmm::create_fmm_matrix(
				p2p_0, p2m_0, p2l_0, m2p_bm, l2p_bm,
				m2m_op, l2l_op, m2l_op,
				tree, lists, std::true_type());


		// integrate operators over fields
		std::cout << "creating integral operators" << std::endl;
		size_t quadrature_order = 10;

		// x2p operators
		fmm::x2p_integral<decltype(m2p_op_0), test_field_t> im2p_op_0(m2p_op_0, quadrature_order);
		fmm::x2p_integral<decltype(l2p_op_0), test_field_t> il2p_op_0(l2p_op_0, quadrature_order);

		// p2p operators
		fmm::p2p_integral<decltype(p2p_op_00), test_field_t, trial_field_t> ip2p_op_00(p2p_op_00, false);
		fmm::p2p_integral<decltype(p2p_op_01), test_field_t, trial_field_t> ip2p_op_01(p2p_op_01, false);

		std::cout << "creating indexed operators" << std::endl;
		auto m2p_0 = fmm::create_x2p_indexed(im2p_op_0, test_space.field_begin<test_field_t>(), test_space.field_end<test_field_t>());
		auto l2p_0 = fmm::create_x2p_indexed(il2p_op_0, test_space.field_begin<test_field_t>(), test_space.field_end<test_field_t>());

		auto p2p_0 = fmm::create_p2x_indexed(
			fmm::create_x2p_indexed(ip2p_op_00,
				test_space.field_begin<test_field_t>(),
				test_space.field_end<test_field_t>()),
			trial_space.field_begin<trial_field_t>(),
			trial_space.field_end<trial_field_t>()
		);
		auto p2p_1 = fmm::create_p2x_indexed(
			fmm::create_x2p_indexed(ip2p_op_01,
				test_space.field_begin<test_field_t>(),
				test_space.field_end<test_field_t>()),
			trial_space.field_begin<trial_field_t>(),
			trial_space.field_end<trial_field_t>()
		);

		std::cout << "startin level data init" << std::endl;
		fmm.init_level_data(tree, 3.0);
		for (size_t c = 0; c < tree.get_n_clusters(); ++c)
			tree[c].set_p_level_data(&fmm.get_level_data(tree[c].get_level()));


		// compute rhs with fmbem
		cvector_t p_field;

		{
			std::cout << "Starting assembling DLP " << std::endl;
			auto dlp_matrix = fmm::create_fmm_matrix(
				p2p_1, p2m_1, p2l_1, m2p_0, l2p_0,
				m2m_op, l2l_op, m2l_op,
				tree, lists, std::false_type());
			std::cout << "DLP assembled" << std::endl;


			// compute solution iteratively
			fmm::matrix_free<decltype(dlp_matrix)> M(dlp_matrix);

			Eigen::GMRES< fmm::matrix_free<decltype(dlp_matrix)>, Eigen::IdentityPreconditioner> solver(M);
			solver.setTolerance(1e-6);
			solver.set_restart(10000);
			p_surf = solver.solve(rhs);

			// compute error
			std::cout << "#iterations: " << solver.iterations() << std::endl;
			std::cout << "#iteration error: " << solver.error() << std::endl;
			std::cout << "timing results: " << std::endl;
			dlp_matrix.get_timer().print(std::cout);


			// export ps
			std::stringstream ss;
			ss << pattern << "_" << freq << "ps.res";
			export_response(ss.str().c_str(), p_surf, k);

		}
#endif

		// read surface solution
		std::stringstream ss;
		ss << pattern << "_" << freq << "ps.res";
		read_excitation(ss.str().c_str(), p_surf, k);

		std::cout << "surface solution read. Wave number: " << k << std::endl;

		{

			typedef NiHu::dirac_field<NiHu::field_view<NiHu::quad_1_elem, NiHu::field_option::constant> > test_field_t;
			auto field_mesh = NiHu::read_off_mesh(field_mesh_name, NiHu::quad_1_tag());
			auto const &test_space = NiHu::dirac(NiHu::constant_view(field_mesh));

			double D = 2.8;
			size_t depth = unsigned(1 + std::log2(k*D));
#ifdef GAUSS
			cluster_tree_t tree(
				fmm::create_field_center_iterator(trial_space.field_begin<trial_field_t>()),
				fmm::create_field_center_iterator(trial_space.field_end<trial_field_t>()),
				fmm::create_elem_center_iterator(field_mesh.begin<elem_t>()),
				fmm::create_elem_center_iterator(field_mesh.end<elem_t>()),
				fmm::divide_depth(depth));
#else
			cluster_tree_t tree(
				fmm::create_elem_center_iterator(mesh.begin<elem_t>()),
				fmm::create_elem_center_iterator(mesh.end<elem_t>()),
				fmm::create_elem_center_iterator(field_mesh.begin<elem_t>()),
				fmm::create_elem_center_iterator(field_mesh.end<elem_t>()),
				fmm::divide_depth(depth));
#endif

			std::cout << tree << std::endl;

			// create interaction lists
			fmm::interaction_lists lists(tree);

			// integrate operators over fields
			std::cout << "creating integral operators" << std::endl;

			// x2p operators
			fmm::x2p_integral<decltype(m2p_op_0), test_field_t> im2p_op_0(m2p_op_0, quadrature_order);
			fmm::x2p_integral<decltype(l2p_op_0), test_field_t> il2p_op_0(l2p_op_0, quadrature_order);

			// p2p operators
			fmm::p2p_integral<decltype(p2p_op_00), test_field_t, trial_field_t> ip2p_op_00(p2p_op_00, false);
			fmm::p2p_integral<decltype(p2p_op_01), test_field_t, trial_field_t> ip2p_op_01(p2p_op_01, false);

			std::cout << "creating indexed operators" << std::endl;
			auto m2p_0 = fmm::create_x2p_indexed(im2p_op_0, test_space.field_begin<test_field_t>(), test_space.field_end<test_field_t>());
			auto l2p_0 = fmm::create_x2p_indexed(il2p_op_0, test_space.field_begin<test_field_t>(), test_space.field_end<test_field_t>());

			auto p2p_0 = fmm::create_p2x_indexed(
				fmm::create_x2p_indexed(ip2p_op_00,
					test_space.field_begin<test_field_t>(),
					test_space.field_end<test_field_t>()),
				trial_space.field_begin<trial_field_t>(),
				trial_space.field_end<trial_field_t>()
			);
			auto p2p_1 = fmm::create_p2x_indexed(
				fmm::create_x2p_indexed(ip2p_op_01,
					test_space.field_begin<test_field_t>(),
					test_space.field_end<test_field_t>()),
				trial_space.field_begin<trial_field_t>(),
				trial_space.field_end<trial_field_t>()
			);

			std::cout << "startin level data init" << std::endl;
			fmm.init_level_data(tree, 3.0);
			for (size_t c = 0; c < tree.get_n_clusters(); ++c)
				tree[c].set_p_level_data(&fmm.get_level_data(tree[c].get_level()));
			std::cout << "level data init finished" << std::endl;


			// compute rhs with fmbem
			cvector_t p_field;

			{
				std::cout << "Starting assembling DLP " << std::endl;
				auto dlp_matrix = fmm::create_fmm_matrix(
					p2p_1, p2m_1, p2l_1, m2p_0, l2p_0,
					m2m_op, l2l_op, m2l_op,
					tree, lists, std::false_type());
				std::cout << "DLP assembled" << std::endl;

				std::cout << "Computing MVP " << std::endl;
				p_field = dlp_matrix * p_surf;
				std::cout << "MVP ready" << std::endl;
			}

			{
				std::cout << "Starting assembling SLP " << std::endl;
				auto slp_matrix = fmm::create_fmm_matrix(
					p2p_0, p2m_0, p2l_0, m2p_0, l2p_0,
					m2m_op, l2l_op, m2l_op,
					tree, lists, std::false_type());
				std::cout << "SLP assembled" << std::endl;
				std::cout << "Computing MVP " << std::endl;
				p_field -= slp_matrix * xct;
				std::cout << "MVP ready" << std::endl;
			}

			std::stringstream ss;
			ss << pattern << "_" << freq << "pf.res";
			export_response(ss.str().c_str(), p_field, k);

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
