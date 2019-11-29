/**
 * \file covariance_2d_bbfmm_matlab.mex.cpp
 * \brief Black box FMM of covariance kernel in 2D with Matlab interface
 * \ingroup tut_fmm
 */

#include "library/covariance_kernel.hpp"

#include "core/field.hpp"
#include "core/function_space.hpp"
#include "fmm/chebyshev_cluster.hpp"
#include "fmm/black_box_fmm.hpp"
#include "fmm/cluster_tree.hpp"
#include "fmm/divide.hpp"
#include "fmm/elem_center_iterator.hpp"
#include "fmm/fmm_indexed.hpp"
#include "fmm/fmm_integrated.hpp"
#include "fmm/fmm_matrix.hpp"
#include "fmm/fmm_precompute.hpp"
#include "library/lib_element.hpp"
#include "util/mex_matrix.hpp"

#include <boost/math/constants/constants.hpp>

#include "mex.h"

#include <cstdlib>
#include <sstream>

typedef NiHu::quad_1_volume_elem elem_t;
typedef NiHu::type2tag<elem_t>::type elem_tag_t;
typedef NiHu::field_view<elem_t, NiHu::field_option::constant> trial_field_t;
typedef NiHu::type2tag<trial_field_t>::type trial_field_tag_t;

typedef trial_field_t test_field_t;
typedef NiHu::type2tag<test_field_t>::type test_field_tag_t;

typedef NiHu::covariance_kernel<NiHu::space_2d<> > kernel_t;
typedef NiHu::fmm::black_box_fmm<kernel_t> fmm_t;

typedef NiHu::mesh<tmp::vector<elem_t> > mesh_t;

typedef fmm_t::cluster_t cluster_t;
typedef NiHu::fmm::cluster_tree<cluster_t> cluster_tree_t;

typedef NiHu::fmm::p2p_precompute<double, 1, 1> p2p_pre_t;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dmatrix_t;
typedef NiHu::fmm::p2x_precompute<dmatrix_t, NiHu::fmm::p2m_tag> p2m_pre_t;
typedef NiHu::fmm::p2x_precompute<dmatrix_t, NiHu::fmm::p2l_tag> p2l_pre_t;
typedef NiHu::fmm::x2p_precompute<dmatrix_t, NiHu::fmm::m2p_tag> m2p_pre_t;
typedef NiHu::fmm::x2p_precompute<dmatrix_t, NiHu::fmm::l2p_tag> l2p_pre_t;
typedef NiHu::fmm::x2x_precompute<fmm_t::m2m::result_t, cluster_t, NiHu::fmm::m2m_tag> m2m_pre_t;
typedef NiHu::fmm::x2x_precompute<fmm_t::l2l::result_t, cluster_t, NiHu::fmm::l2l_tag> l2l_pre_t;
typedef NiHu::fmm::x2x_precompute<fmm_t::m2l::result_t, cluster_t, NiHu::fmm::m2l_tag> m2l_pre_t;

typedef NiHu::fmm::fmm_matrix<
	p2p_pre_t,
	p2m_pre_t,
	p2l_pre_t,
	m2p_pre_t,
	l2p_pre_t,
	m2m_pre_t,
	l2l_pre_t,
	m2l_pre_t> fmm_matrix_t;

// Mex matrix types
typedef NiHu::mex::real_matrix<double> dMatrix;


class fmm_matlab
{
public:
	fmm_matlab()
		: p_surf_mesh(nullptr)
		, p_tree(nullptr)
		, p_lists(nullptr)
		, p_fmm(nullptr)
		, p_fmm_matrix(nullptr)
		, p_I_pre(nullptr)
	{
	}

	void create_mesh(dMatrix const &surf_nodes, dMatrix const &surf_elems)
	{
		p_surf_mesh = new mesh_t(NiHu::create_mesh(surf_nodes, surf_elems, NiHu::quad_1_volume_tag()));
	}

	template <class DivideDerived>
	void create_tree(NiHu::fmm::divide_base<DivideDerived> const &divide)
	{
		p_tree = new cluster_tree_t(
			NiHu::fmm::create_elem_center_iterator(p_surf_mesh->begin<elem_t>()),
			NiHu::fmm::create_elem_center_iterator(p_surf_mesh->end<elem_t>()),
			divide
		);

		// create interaction lists
		p_lists = new NiHu::fmm::interaction_lists(*p_tree);
	}

	void create_matrix()
	{

		kernel_t kernel(m_sigma, m_cov_length);
		p_fmm = new fmm_t(kernel);

		// initialize tree data
		for (size_t c = 0; c < p_tree->get_n_clusters(); ++c)
			(*p_tree)[c].set_chebyshev_order(m_cheb_order);

		size_t far_field_quadrature_order = 5;

		// create functors
		auto int_fctr = NiHu::fmm::create_integrated_functor(test_field_tag_t(), trial_field_tag_t(),
			far_field_quadrature_order, true);

		auto const &trial_space = NiHu::constant_view(*p_surf_mesh);
		auto const &test_space = trial_space;

		auto idx_fctr = create_indexed_functor(
			test_space.template field_begin<test_field_t>(),
			test_space.template field_end<test_field_t>(),
			trial_space.template field_begin<trial_field_t>(),
			trial_space.template field_end<trial_field_t>(),
			*p_tree);

		auto pre_fctr = NiHu::fmm::create_precompute_functor(*p_tree, *p_lists);

		p_fmm_matrix = new fmm_matrix_t(NiHu::fmm::create_fmm_matrix(
			pre_fctr(idx_fctr(int_fctr(p_fmm->create_p2p()))),
			pre_fctr(idx_fctr(int_fctr(p_fmm->create_p2m()))),
			pre_fctr(idx_fctr(int_fctr(p_fmm->create_p2l()))),
			pre_fctr(idx_fctr(int_fctr(p_fmm->create_m2p()))),
			pre_fctr(idx_fctr(int_fctr(p_fmm->create_l2p()))),
			pre_fctr(idx_fctr(p_fmm->create_m2m())),
			pre_fctr(idx_fctr(p_fmm->create_l2l())),
			pre_fctr(idx_fctr(p_fmm->create_m2l())),
			*p_tree,
			*p_lists));

		p_I_pre = new p2p_pre_t(
			pre_fctr(
				idx_fctr(
					NiHu::fmm::create_identity_p2p_integral(
						test_field_tag_t(),
						trial_field_tag_t()
					)
				)
			)
			);
	}

	template <class LhsDerived, class RhsDerived>
	void mvp(Eigen::MatrixBase<LhsDerived> &&res, Eigen::MatrixBase<RhsDerived> const &src)
	{
		res = (*p_fmm_matrix) * src;
	}

	void print_tree() const
	{
		// debug output
		std::stringstream ss;
		ss << *p_tree;

		mexPrintf("Tree:\n%s\n", ss.str().c_str());
	}

	~fmm_matlab()
	{
		delete p_I_pre;
		delete p_fmm_matrix;
		delete p_fmm;
		delete p_lists;
		delete p_tree;
		delete p_surf_mesh;
	}

	void set_sigma(double sigma)
	{
		m_sigma = sigma;
	}

	void set_cov_length(double length)
	{
		m_cov_length = length;
	}

	void set_cheb_order(size_t order)
	{
		m_cheb_order = order;
	}

	Eigen::SparseMatrix<double> const &get_sparse_identity() const
	{
		return p_I_pre->get_sparse_matrix();
	}

private:
	mesh_t *p_surf_mesh;
	cluster_tree_t *p_tree;
	NiHu::fmm::interaction_lists *p_lists;
	fmm_t *p_fmm;
	fmm_matrix_t *p_fmm_matrix;
	p2p_pre_t *p_I_pre;

	double m_sigma;
	double m_cov_length;
	size_t m_cheb_order;
};

fmm_matlab *p = nullptr;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
	char *input_option;

	/* input must be a string */
	if (mxIsChar(prhs[0]) != 1)
		mexErrMsgIdAndTxt("MATLAB:revord:inputNotString",
			"Input must be a string.");

	input_option = mxArrayToString(prhs[0]);

	if (!strcmp(input_option, "init"))
	{
		if (p != nullptr)
		{
			mexErrMsgIdAndTxt("NiHu:covariance_2d_bbfmm_matlab:invalid_state",
				"Command \"%s\" called in invalid state", input_option);
		}
		else
		{
			p = new fmm_matlab;
		}
	}

	else if (!strcmp(input_option, "set"))
	{
		int n_pairs = (nrhs - 1) / 2;
		for (int i = 0; i < n_pairs; ++i)
		{
			if (mxIsChar(prhs[2 * i + 1]) != 1)
			{
				mexErrMsgIdAndTxt("NiHu:covariance_2d_bbfmm_matlab:invalid_input",
					"Parameter name must be a string for the command \"%s\".", input_option);
			}
			char const *what_to_set = mxArrayToString(prhs[2 * i + 1]);

			if (!strcmp(what_to_set, "sigma"))
			{
				double sigma = mxGetScalar(prhs[2 * i + 2]);
				p->set_sigma(sigma);
			}
			else if (!strcmp(what_to_set, "cov_length"))
			{
				double cov_length = mxGetScalar(prhs[2 * i + 2]);
				p->set_cov_length(cov_length);
			}
			else if (!strcmp(what_to_set, "cheb_order"))
			{
				double cheb_order = mxGetScalar(prhs[2 * i + 2]);
				p->set_cheb_order(size_t(cheb_order));
			}
			else
			{
				mexErrMsgIdAndTxt("NiHu:covariance_2d_bbfmm_matlab:invalid_parameter",
					"Unknown input parameter: \"%s\"", what_to_set);
			}
		}
	}

	else if (!strcmp(input_option, "mesh"))
	{
		if (p == nullptr)
		{
			mexErrMsgIdAndTxt("NiHu:covariance_2d_bbfmm_matlab:invalid_state",
				"Command \"%s\" called in invalid state", input_option);
		}
		else
		{
			p->create_mesh(dMatrix(prhs[1]), dMatrix(prhs[2]));
		}
	}

	// "tree" command - Build the tree
	else if (!strcmp(input_option, "tree"))
	{
		if (p == nullptr)
		{
			mexErrMsgIdAndTxt("NiHu:covariance_2d_bbfmm_matlab:invalid_state",
				"Command \"%s\" called in invalid state", input_option);
		}
		else
		{
			if (mxIsChar(prhs[1]) != 1)
			{
				mexErrMsgIdAndTxt("NiHu:covariance_2d_bbfmm_matlab:invalid_input",
					"Tree division method must be a string for the command \"%s\".", input_option);
			}
			char const *divide_option = mxArrayToString(prhs[1]);
			if (!strcmp(divide_option, "divide_depth"))
			{
				p->create_tree(NiHu::fmm::divide_depth(size_t(mxGetScalar(prhs[2]))));
			}
			else if (!strcmp(divide_option, "divide_num_nodes"))
			{
				p->create_tree(NiHu::fmm::divide_num_nodes(size_t(mxGetScalar(prhs[2]))));
			}
			else if (!strcmp(divide_option, "divide_diameter"))
			{
				p->create_tree(NiHu::fmm::divide_diameter(mxGetScalar(prhs[2])));
			}
			else
			{
				mexErrMsgIdAndTxt("NiHu:covariance_2d_bbfmm_matlab:invalid_divide_option",
					"Unknown divide option: \"%s\"", divide_option);
			}
		}
	}

	// "matrix" option - assemble FMM matrix
	else if (!strcmp(input_option, "matrix"))
	{
		if (p == nullptr)
		{
			mexErrMsgIdAndTxt("NiHu:covariance_2d_bbfmm_matlab:invalid_state",
				"Command \"%s\" called in invalid state", input_option);
		}
		else
		{
			p->create_matrix();
		}
	}

	else if (!strcmp(input_option, "get_sparse_identity"))
	{
	if (p == nullptr)
	{
		mexErrMsgIdAndTxt("NiHu:covariance_2d_bbfmm_matlab:invalid_state",
			"Command \"%s\" called in invalid state", input_option);
	}
	else
	{
		Eigen::SparseMatrix<double> const &mat = p->get_sparse_identity();
		plhs[0] = mxCreateSparse(mat.rows(), mat.cols(), mwSize(mat.nonZeros()), mxREAL);
		
		mwIndex *ridx = mxGetIr(plhs[0]);
		mwIndex *cidx = mxGetJc(plhs[0]);
		int c = 0;
		int k = 0;
		double *v = mxGetPr(plhs[0]);
		for (k = 0; k < mat.outerSize(); ++k) {
			cidx[k] = c;
			for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it) {
				v[c] = it.value();
				ridx[c] = it.row();
				++c;
			}
		}
		cidx[k] = c;
	}
 }

	else if (!strcmp(input_option, "mvp"))
	{
		if (p == nullptr)
		{
			mexErrMsgIdAndTxt("NiHu:covariance_2d_bbfmm_matlab:invalid_state",
				"Command \"%s\" called in invalid state", input_option);
		}
		else
		{
			dMatrix xct(prhs[1]);
			dMatrix res(xct.rows(), 1, plhs[0]);
			p->mvp(res.col(0), xct.col(0));
		}
	}

	// "cleanup" - destroy the main object
	else if (!strcmp(input_option, "cleanup"))
	{
		if (p == nullptr)
		{
			mexErrMsgIdAndTxt("NiHu:covariance_2d_bbfmm_matlab:invalid_state",
				"Command \"%s\" called in invalid state", input_option);
		}
		else
		{
			delete p;
			p = nullptr;
		}
	}

	// The option was not valid
	else
	{
		mexErrMsgIdAndTxt("NiHu:covariance_2d_bbfmm_matlab:invalid_option",
			"Unknown input option: \"%s\"", input_option);
	}
}
