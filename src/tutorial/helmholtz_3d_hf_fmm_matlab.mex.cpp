#include "core/field.hpp"
#include "core/function_space.hpp"
#include "fmm/divide.hpp"
#include "fmm/elem_center_iterator.hpp"
#include "fmm/fmm_integrated.hpp"
#include "fmm/fmm_indexed.hpp"
#include "fmm/fmm_matrix.hpp"
#include "fmm/fmm_operator_collection.hpp"
#include "fmm/fmm_precompute.hpp"
#include "fmm/helmholtz_3d_hf_fmm.hpp"

#include "library/lib_element.hpp"
#include "util/mex_matrix.hpp"

#include <boost/math/constants/constants.hpp>

#include "mex.h"

#include <cstdlib>
#include <sstream>

typedef NiHu::quad_1_elem elem_t;
typedef NiHu::type2tag<elem_t>::type elem_tag_t;
typedef NiHu::field_view<elem_t, NiHu::field_option::constant> trial_field_t;
typedef NiHu::type2tag<trial_field_t>::type trial_field_tag_t;

typedef NiHu::dirac_field<trial_field_t> test_field_t;
typedef NiHu::type2tag<test_field_t>::type test_field_tag_t;

typedef NiHu::fmm::helmholtz_3d_hf_fmm<double> fmm_t;

typedef NiHu::mesh<tmp::vector<elem_t> > mesh_t;

typedef fmm_t::cluster_t cluster_t;
typedef NiHu::fmm::cluster_tree<cluster_t> cluster_tree_t;

typedef NiHu::fmm::p2p_precompute<std::complex<double>, 1, 1> p2p_t;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dmatrix_t;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cmatrix_t;

typedef NiHu::fmm::p2x_precompute<cmatrix_t, NiHu::fmm::p2m_tag> p2m_t;
typedef NiHu::fmm::p2x_precompute<cmatrix_t, NiHu::fmm::p2l_tag> p2l_t;
typedef NiHu::fmm::x2p_precompute<cmatrix_t, NiHu::fmm::m2p_tag> m2p_t;
typedef NiHu::fmm::x2p_precompute<cmatrix_t, NiHu::fmm::l2p_tag> l2p_t;
typedef NiHu::fmm::x2x_precompute<fmm_t::m2m::result_t, cluster_t, NiHu::fmm::m2m_tag> m2m_t;
typedef NiHu::fmm::x2x_precompute<fmm_t::l2l::result_t, cluster_t, NiHu::fmm::l2l_tag> l2l_t;
typedef NiHu::fmm::x2x_precompute<fmm_t::m2l::result_t, cluster_t, NiHu::fmm::m2l_tag> m2l_t;

typedef NiHu::fmm::fmm_matrix<
	p2p_t, p2m_t, p2l_t,
	m2p_t, l2p_t, 
	m2m_t, l2l_t, m2l_t> fmm_matrix_t;



// Mex matrix types
typedef NiHu::mex::real_matrix<double> dmex_matrix_t;
typedef NiHu::mex::complex_matrix<double> cmex_matrix_t;


class fmm_matlab
{
public:
	fmm_matlab() 
		: p_surf_mesh(nullptr)
		, p_field_mesh(nullptr)
		, p_tree(nullptr)
		, p_lists(nullptr)
		, p_fmm(nullptr)
		, p_slp_matrix(nullptr)
		, p_dlp_matrix(nullptr)
	{
		
	}
	
	void create_mesh(dmex_matrix_t const& surf_nodes, dmex_matrix_t const& surf_elems)
	{
		p_surf_mesh = new mesh_t(NiHu::create_mesh(surf_nodes, surf_elems, NiHu::quad_1_tag()));
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
	
	void create_matrices()
	{
		// create the FMM instance
		p_fmm = new fmm_t(m_wave_number);
		p_fmm->set_accuracy(m_accuracy);

		// initialize tree data
		std::cout << "Initializing tree data ..." << std::endl;
		p_fmm->init_level_data(*p_tree);
		for (size_t c = 0; c < p_tree->get_n_clusters(); ++c)
			(*p_tree)[c].set_p_level_data(&p_fmm->get_level_data((*p_tree)[c].get_level()));

		
		size_t far_field_quadrature_order = 5;

		// create functors
		auto int_fctr = NiHu::fmm::create_integrated_functor(
			test_field_tag_t(), trial_field_tag_t(), far_field_quadrature_order, true);

		auto const &trial_space = NiHu::constant_view(*p_surf_mesh);
		auto const &test_space  = NiHu::dirac(trial_space);

		auto idx_fctr = NiHu::fmm::create_indexed_functor(
			test_space.template field_begin<test_field_t>(),
			test_space.template field_end<test_field_t>(),
			trial_space.template field_begin<trial_field_t>(),
			trial_space.template field_end<trial_field_t>(),
			*p_tree);

		auto pre_fctr = NiHu::fmm::create_precompute_functor(*p_tree, *p_lists);

		// Burton-Miller coupling constant
		std::complex<double> alpha(0.0, -1.0 / m_wave_number);

		// integration
		auto I = NiHu::fmm::create_identity_p2p_integral(test_field_tag_t(), trial_field_tag_t());
		// Create the operator collection for the LHS matrix
		auto dlp_collection = NiHu::fmm::create_fmm_operator_collection(
			int_fctr(p_fmm->template create_p2p<0, 1>())
				+ alpha * int_fctr(p_fmm->template create_p2p<1, 1>())
				- 0.5 * I,
			int_fctr(p_fmm->template create_p2m<1>()),
			int_fctr(p_fmm->template create_p2l<1>()),
			int_fctr(p_fmm->template create_m2p<0>())
				+ alpha * int_fctr(p_fmm->template create_m2p<1>()),
			int_fctr(p_fmm->template create_l2p<0>())
				+ alpha * int_fctr(p_fmm->template create_l2p<1>()),
			p_fmm->create_m2m(),
			p_fmm->create_m2l(),
			p_fmm->create_l2l()
		);

		// Create the collection of additional operators for the RHS matrix
		auto slp_collection = NiHu::fmm::create_fmm_operator_collection(
			int_fctr(p_fmm->template create_p2p<0, 0>())
				+ alpha * int_fctr(p_fmm->template create_p2p<1, 0>())
				+ (alpha / 2.0) * I,
			int_fctr(p_fmm->template create_p2m<0>()),
			int_fctr(p_fmm->template create_p2l<0>())
		);

		// indexing
		auto dlp_cix_collection = dlp_collection.transform(idx_fctr);
		auto slp_cix_collection = slp_collection.transform(idx_fctr);

		// precomputation
		std::cout << "Precomputing fmm operators ..." << std::endl;
		
		auto dlp_pre_collection = dlp_cix_collection.transform(pre_fctr);
		auto slp_pre_collection = slp_cix_collection.transform(pre_fctr);

		// create slp matrix object
		std::cout << "Assembling slp matrix ..." << std::endl;
		p_slp_matrix = new fmm_matrix_t(NiHu::fmm::create_fmm_matrix(
			slp_pre_collection.get(NiHu::fmm::p2p_tag()),
			slp_pre_collection.get(NiHu::fmm::p2m_tag()),
			slp_pre_collection.get(NiHu::fmm::p2l_tag()),
			dlp_pre_collection.get(NiHu::fmm::m2p_tag()),
			dlp_pre_collection.get(NiHu::fmm::l2p_tag()),
			dlp_pre_collection.get(NiHu::fmm::m2m_tag()),
			dlp_pre_collection.get(NiHu::fmm::l2l_tag()),
			dlp_pre_collection.get(NiHu::fmm::m2l_tag()),
			*p_tree, *p_lists));

		// create matrix object
		std::cout << "Assembling dlp matrix ..." << std::endl;
		p_dlp_matrix = new fmm_matrix_t(create_fmm_matrix(
			dlp_pre_collection,
			*p_tree, *p_lists));
	}
	
	template <class LhsDerived, class RhsDerived>
	void mvp_dlp(Eigen::MatrixBase<LhsDerived> &&res, Eigen::MatrixBase<RhsDerived> const &src)
	{
		res = (*p_dlp_matrix) * src;
	}
	
	template <class LhsDerived, class RhsDerived>
	void mvp_slp(Eigen::MatrixBase<LhsDerived> &&res, Eigen::MatrixBase<RhsDerived> const &src)
	{
		res = (*p_slp_matrix) * src;
	}
	
	void print_tree()
	{
		if (!p_tree)
			return;
		
		// debug output
		std::stringstream ss;
		ss << *p_tree;
		
		mexPrintf("Tree:\n%s\n", ss.str().c_str());
	}
	
	void set_accuracy(double accuracy)
	{
		m_accuracy = accuracy;
	}
	
	void set_wave_number(double wave_number)
	{
		m_wave_number = wave_number;
	}
	
	
	~fmm_matlab()
	{
		delete p_surf_mesh;
		delete p_field_mesh;
		delete p_tree;
		delete p_lists;
		delete p_fmm;
		delete p_slp_matrix;
		delete p_dlp_matrix;
		
	}
private:
	mesh_t *p_surf_mesh;
	mesh_t *p_field_mesh;
	cluster_tree_t *p_tree;
	NiHu::fmm::interaction_lists *p_lists;
	fmm_t *p_fmm;
	fmm_matrix_t *p_slp_matrix;
	fmm_matrix_t *p_dlp_matrix;
	
	double m_wave_number;
	double m_accuracy;
	
};

fmm_matlab *p = nullptr;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
	char const *input_option;
	
	 /* input must be a string */
    if ( mxIsChar(prhs[0]) != 1) {
		mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_surf:invalid_input",
              "First input parameter must be a string.");
	}
	input_option = mxArrayToString(prhs[0]);
	
	if (!strcmp(input_option, "init")) 
	{
		if (p != nullptr) {
			mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_matlab:invalid_state",
				"Command \"init\" called in invalid state");
		} else {
			p = new fmm_matlab;
		}
	}
	
	else if (!strcmp(input_option, "set"))
	{
		int n_pairs = (nrhs - 1) / 2;
		// Go through parameter pairs
		for (int i = 0; i < n_pairs; ++i)
		{
			if ( mxIsChar(prhs[2*i + 1]) != 1) {
				mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_surf:invalid_input",
				"Parameter name must be a string for the command \"set\".");
			}
			
			char const *what_to_set = mxArrayToString(prhs[2 * i + 1]);
			if (!strcmp(what_to_set, "wave_number"))
			{
				double k = mxGetScalar(prhs[2 * i + 2]);
				p->set_wave_number(k);
			}
			if (!strcmp(what_to_set, "accuracy"))
			{
				double accuracy = mxGetScalar(prhs[2 * i + 2]);
				p->set_accuracy(accuracy);
			}
			else
			{
				mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_surf:invalid_parameter",
					"Unknown parameter name to set: \"%s\"", what_to_set);
			}
		}
	}
	
	else if (!strcmp(input_option, "mesh"))
	{
		if (p == nullptr)
		{
			mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_matlab:invalid_state",
				"Command \"mesh\" called in invalid state");
		}
		else
		{
			p->create_mesh(dmex_matrix_t(prhs[1]), dmex_matrix_t(prhs[2]));
		}
	}
	
	else if (!strcmp(input_option, "tree"))
	{
		if (p == nullptr)
		{
			mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_matlab:invalid_state",
				"Command \"tree\" called in invalid state");
		}
		else
		{
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
	
	else if (!strcmp(input_option, "matrix"))
	{
		if (p == nullptr)
		{
			mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_matlab:invalid_state",
				"Command \"matrix\" called in invalid state");
		}
		else
		{
			p->create_matrices();
		}
	}
	
	else if (!strcmp(input_option, "mvp_dlp"))
	{
		if (p == nullptr)
		{
			mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_matlab:invalid_state",
				"Command \"%s\" called in invalid state", input_option);
		}
		else
		{
			cmex_matrix_t xct(prhs[1]);
			cmex_matrix_t res(xct.rows(), 1, plhs[0]);
			p->mvp_dlp(res.col(0), xct.col(0));
		}
	}
	
	else if (!strcmp(input_option, "mvp_slp"))
	{
		if (p == nullptr)
		{
			mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_matlab:invalid_state",
				"Command \"%s\" called in invalid state", input_option);
		}
		else
		{
			cmex_matrix_t xct(prhs[1]);
			cmex_matrix_t res(xct.rows(), 1, plhs[0]);
			p->mvp_slp(res.col(0), xct.col(0));
		}
	}
	
	
	
	// Cleanup option 
	else if (!strcmp(input_option, "cleanup"))
	{
		if (p == nullptr)
		{
			mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_matlab:invalid_state",
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
		mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_matlab:invalid_option",
			"Unknown input option: \"%s\"", input_option);
	}	

}