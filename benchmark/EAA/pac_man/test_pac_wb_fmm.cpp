#include "cluster_tree.hpp"
#include "divide.hpp"
#include "elem_center_iterator.hpp"
#include "fmm_matrix.hpp"
#include "helmholtz_2d_wb_fmm.hpp"
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

#include <Eigen/IterativeLinearSolvers>
#include "GMRES.h"

// basic type parameter inputs
typedef NiHu::line_1_tag tag_t;
typedef double wave_number_t;

// computing the fmm type
typedef fmm::helmholtz_2d_wb_fmm<wave_number_t> fmm_t;

// computing the fmbem type
typedef NiHu::tag2element<tag_t>::type elem_t;
typedef NiHu::field_view<elem_t, NiHu::field_option::constant> trial_field_t;
typedef NiHu::dirac_field<trial_field_t> test_field_t;
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

void export_response(std::string fname, cvector_t const &res, double k)
{
	std::ofstream ofs(fname);
	ofs << k << '\n';
	ofs << res.rows() << '\n';
	for (size_t i = 0; i < res.rows(); ++i)
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

	// read mesh
	auto surf_mesh = NiHu::read_off_mesh(surf_mesh_name, tag_t());

	// read field
	auto field_mesh = NiHu::read_off_mesh(field_mesh_name, tag_t());

	// read excitation
	double k;
	cvector_t xct, p_surf;
	read_excitation(xct_name, xct, k);

	std::cout << "wave number: " << k << std::endl;
	// create kernel, fmm and fmbem objects
	fmm_t fmm(k);

	// solve surface system
	{
		auto const &trial_space = NiHu::constant_view(surf_mesh);
		auto const &test_space = NiHu::dirac(trial_space);

		cluster_tree_t tree(
			fmm::create_elem_center_iterator(surf_mesh.begin<elem_t>()),
			fmm::create_elem_center_iterator(surf_mesh.end<elem_t>()),
			fmm::create_elem_center_iterator(surf_mesh.begin<elem_t>()),
			fmm::create_elem_center_iterator(surf_mesh.end<elem_t>()),
			fmm::divide_num_nodes(10));

		std::cout << tree << std::endl;

		// create interaction lists
		fmm::interaction_lists lists(tree);

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

		// p2x operators
		fmm::p2x_integral<decltype(p2m_op_0), trial_field_t> ip2m_op_0(p2m_op_0, quadrature_order);
		fmm::p2x_integral<decltype(p2l_op_0), trial_field_t> ip2l_op_0(p2l_op_0, quadrature_order);
		fmm::p2x_integral<decltype(p2m_op_1), trial_field_t> ip2m_op_1(p2m_op_1, quadrature_order);
		fmm::p2x_integral<decltype(p2l_op_1), trial_field_t> ip2l_op_1(p2l_op_1, quadrature_order);

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

		auto p2m_0 = fmm::create_p2x_indexed(ip2m_op_0, trial_space.field_begin<trial_field_t>(), trial_space.field_end<trial_field_t>());
		auto p2l_0 = fmm::create_p2x_indexed(ip2l_op_0, trial_space.field_begin<trial_field_t>(), trial_space.field_end<trial_field_t>());
		auto p2m_1 = fmm::create_p2x_indexed(ip2m_op_1, trial_space.field_begin<trial_field_t>(), trial_space.field_end<trial_field_t>());
		auto p2l_1 = fmm::create_p2x_indexed(ip2l_op_1, trial_space.field_begin<trial_field_t>(), trial_space.field_end<trial_field_t>());

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

		fmm.init_level_data(tree);
		fmm.print_level_data();
		for (size_t c = 0; c < tree.get_n_clusters(); ++c)
			tree[c].set_p_level_data(&fmm.get_level_data(tree[c].get_level()));

		std::cout << "Starting assembling matrices " << std::endl;

		// create matrix objects
		auto slp_matrix = fmm::create_fmm_matrix(
			p2p_0, p2m_0, p2l_0, m2p_bm, l2p_bm,
			m2m_op, l2l_op, m2l_op,
			tree, lists, std::true_type());

		auto dlp_matrix = fmm::create_fmm_matrix(
			p2p_1, p2m_1, p2l_1, m2p_bm, l2p_bm,
			m2m_op, l2l_op, m2l_op,
			tree, lists, std::true_type());

		std::cout << "Matrices ready " << std::endl;

		// compute rhs with fmbem
		decltype(slp_matrix)::response_t rhs = (slp_matrix * xct + alpha / 2. * xct).eval();

		std::cout << "rhs ready" << std::endl;

		std::cout << "Starting iterative solution" << std::endl;

		// compute solution iteratively
		fmm::matrix_free<decltype(dlp_matrix)> M(dlp_matrix);
		Eigen::GMRES<fmm::matrix_free<decltype(dlp_matrix)>, Eigen::IdentityPreconditioner> solver(M);
		solver.setTolerance(1e-8);
		p_surf = solver.solve(rhs);

		// compute error
		std::cout << "#iterations: " << solver.iterations() << std::endl;
		std::cout << "#iteration error: " << solver.error() << std::endl;
		std::cout << "timing results: " << std::endl;
		dlp_matrix.get_timer().print(std::cout);

		export_response(surf_res_name, p_surf, k);

	}

	// evaluate field pressure
	{
		auto const &trial_space = NiHu::constant_view(surf_mesh);
		auto const &test_space = NiHu::dirac(NiHu::constant_view(field_mesh));

		cluster_tree_t tree(
			fmm::create_elem_center_iterator(surf_mesh.begin<elem_t>()),
			fmm::create_elem_center_iterator(surf_mesh.end<elem_t>()),
			fmm::create_elem_center_iterator(field_mesh.begin<elem_t>()),
			fmm::create_elem_center_iterator(field_mesh.end<elem_t>()),
			fmm::divide_num_nodes(10));

		std::cout << tree << std::endl;

		// create interaction lists
		fmm::interaction_lists lists(tree);

		// get x2x fmm operators
		auto m2m_op = fmm.create_m2m();
		auto l2l_op = fmm.create_l2l();
		auto m2l_op = fmm.create_m2l();

		// p2x operators and derivatives
		auto p2m_op_0 = fmm.create_p2m<0>();
		auto p2l_op_0 = fmm.create_p2l<0>();
		auto p2m_op_1 = fmm.create_p2m<1>();
		auto p2l_op_1 = fmm.create_p2l<1>();

		// x2p operators
		auto m2p_op_0 = fmm.create_m2p<0>();
		auto l2p_op_0 = fmm.create_l2p<0>();

		// p2p operators and derivatives
		auto p2p_op_00 = fmm.create_p2p<0, 0>();
		auto p2p_op_01 = fmm.create_p2p<0, 1>();

		// integrate operators over fields
		std::cout << "creating integral operators" << std::endl;
		size_t quadrature_order = 10;

		// p2x operators
		fmm::p2x_integral<decltype(p2m_op_0), trial_field_t> ip2m_op_0(p2m_op_0, quadrature_order);
		fmm::p2x_integral<decltype(p2l_op_0), trial_field_t> ip2l_op_0(p2l_op_0, quadrature_order);
		fmm::p2x_integral<decltype(p2m_op_1), trial_field_t> ip2m_op_1(p2m_op_1, quadrature_order);
		fmm::p2x_integral<decltype(p2l_op_1), trial_field_t> ip2l_op_1(p2l_op_1, quadrature_order);

		// x2p operators
		fmm::x2p_integral<decltype(m2p_op_0), test_field_t> im2p_op_0(m2p_op_0, quadrature_order);
		fmm::x2p_integral<decltype(l2p_op_0), test_field_t> il2p_op_0(l2p_op_0, quadrature_order);

		// p2p operators
		fmm::p2p_integral<decltype(p2p_op_00), test_field_t, trial_field_t> ip2p_op_00(p2p_op_00, false);
		fmm::p2p_integral<decltype(p2p_op_01), test_field_t, trial_field_t> ip2p_op_01(p2p_op_01, false);

		// create indexed fmbem operators

		std::cout << "creating indexed operators" << std::endl;
		auto p2m_0 = fmm::create_p2x_indexed(ip2m_op_0, trial_space.field_begin<trial_field_t>(), trial_space.field_end<trial_field_t>());
		auto p2l_0 = fmm::create_p2x_indexed(ip2l_op_0, trial_space.field_begin<trial_field_t>(), trial_space.field_end<trial_field_t>());
		auto p2m_1 = fmm::create_p2x_indexed(ip2m_op_1, trial_space.field_begin<trial_field_t>(), trial_space.field_end<trial_field_t>());
		auto p2l_1 = fmm::create_p2x_indexed(ip2l_op_1, trial_space.field_begin<trial_field_t>(), trial_space.field_end<trial_field_t>());

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
		fmm.init_level_data(tree);
		fmm.print_level_data();
		for (size_t c = 0; c < tree.get_n_clusters(); ++c)
			tree[c].set_p_level_data(&fmm.get_level_data(tree[c].get_level()));


		// create matrix objects


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

			slp_matrix.get_timer().print(std::cout);
		}

		export_response(field_res_name, p_field, k);
	}

	return 0;
}
