/**
 * @file helmholtz_3d_hffmm_mex.mex.cpp 
 * @brief Helmholtz fast multipole solver in 3D with MEX interface
 * @ingroup app_helmholtz
 * 
 * @details
 * High frequency fast multipole method for the Helmholtz equation, 
 * collocational formalism, Burton-Miller approach for mitigating fictitious 
 * eigenfrequencies. 
 */

#include "nihu/core/field.hpp"
#include "nihu/core/function_space.hpp"
#include "nihu/fmm/divide.hpp"
#include "nihu/fmm/elem_center_iterator.hpp"
#include "nihu/fmm/fmm_integrated.hpp"
#include "nihu/fmm/fmm_indexed.hpp"
#include "nihu/fmm/fmm_matrix.hpp"
#include "nihu/fmm/fmm_operator_collection.hpp"
#include "nihu/fmm/fmm_precompute.hpp"
#include "nihu/fmm/fmm_timer.h"
#include "nihu/fmm/helmholtz_3d_hf_fmm.hpp"
#include "nihu/library/lib_element.hpp"
#include "nihu/util/mex_matrix.hpp"

#include <boost/math/constants/constants.hpp>

#include "mex.h"

#include <cstdlib>
#include <sstream>

#define NIHU_THIS_MEX_NAME "helmholtz_3d_hffmm_mex"

// Eigen matrix types
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dmatrix_t;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cmatrix_t;

// Mex matrix types
typedef NiHu::mex::real_matrix<double> dmex_matrix_t;
typedef NiHu::mex::complex_matrix<double> cmex_matrix_t;

/**
 * @brief Helper class for storing FMM assembly times 
 */
class fmm_assembly_times
{
	using timer = NiHu::fmm::fmm_timer;
public:
	/** 
	 * @brief Fill the assembly time from an operator collection
	 * @tparam Collection The operator collection class
	 * @param[in] coll Collection instance
	 */
	template <class Collection>
	void fill_times(Collection const & coll)
	{
		m_times[timer::M2M] = coll.get(NiHu::fmm::op_tags::idx2tag<timer::M2M>::type()).get_assembly_time();
		m_times[timer::L2L] = coll.get(NiHu::fmm::op_tags::idx2tag<timer::L2L>::type()).get_assembly_time();
		m_times[timer::M2L] = coll.get(NiHu::fmm::op_tags::idx2tag<timer::M2L>::type()).get_assembly_time();
		
		m_times[timer::P2M] = coll.get(NiHu::fmm::op_tags::idx2tag<timer::P2M>::type()).get_assembly_time();
		m_times[timer::P2L] = coll.get(NiHu::fmm::op_tags::idx2tag<timer::P2L>::type()).get_assembly_time();
		m_times[timer::L2P] = coll.get(NiHu::fmm::op_tags::idx2tag<timer::L2P>::type()).get_assembly_time();
		m_times[timer::M2P] = coll.get(NiHu::fmm::op_tags::idx2tag<timer::M2P>::type()).get_assembly_time();
		
		m_times[timer::P2P] = coll.get(NiHu::fmm::op_tags::idx2tag<timer::P2P>::type()).get_assembly_time();
	}
	
	/**
	 * @brief Returns the assembly time for the given operation index
	 * @return Assembly time in microsecond units
	 */
	size_t const get_time(size_t idx) const
	{
		return m_times[idx];
	}
	
private:
	size_t m_times[timer::NUM_TIME_INDICES];
};

/**
 * @brief Class for FMM-related data and methods 
 */
class fmm_matlab
{
    typedef NiHu::quad_1_elem elem_t;
    typedef NiHu::type2tag<elem_t>::type elem_tag_t;
    typedef NiHu::field_view<elem_t, NiHu::field_option::constant> trial_field_t;
    typedef NiHu::type2tag<trial_field_t>::type trial_field_tag_t;

    typedef NiHu::dirac_field<trial_field_t> test_field_t;
    typedef NiHu::type2tag<test_field_t>::type test_field_tag_t;

    typedef double wavenumber_t;
    typedef NiHu::fmm::helmholtz_3d_hf_fmm<wavenumber_t> fmm_t;

    typedef NiHu::mesh<tmp::vector<elem_t> > mesh_t;

    typedef fmm_t::cluster_t cluster_t;
    typedef NiHu::fmm::cluster_tree<cluster_t> cluster_tree_t;

    typedef NiHu::fmm::p2p_precompute<std::complex<double>, 1, 1> p2p_t;

    // Typedefs for the precomputed operators
    typedef NiHu::fmm::p2x_precompute<cmatrix_t, NiHu::fmm::p2m_tag> p2m_t;
    typedef NiHu::fmm::p2x_precompute<cmatrix_t, NiHu::fmm::p2l_tag> p2l_t;
    typedef NiHu::fmm::x2p_precompute<cmatrix_t, NiHu::fmm::m2p_tag> m2p_t;
    typedef NiHu::fmm::x2p_precompute<cmatrix_t, NiHu::fmm::l2p_tag> l2p_t;
    typedef NiHu::fmm::x2x_precompute<fmm_t::m2m::result_t, cluster_t, NiHu::fmm::m2m_tag> m2m_t;
    typedef NiHu::fmm::x2x_precompute<fmm_t::l2l::result_t, cluster_t, NiHu::fmm::l2l_tag> l2l_t;
    typedef NiHu::fmm::x2x_precompute<fmm_t::m2l::result_t, cluster_t, NiHu::fmm::m2l_tag> m2l_t;

    // The FMM matrix type
    typedef NiHu::fmm::fmm_matrix<
        p2p_t, p2m_t, p2l_t,
        m2p_t, l2p_t, 
        m2m_t, l2l_t, m2l_t> fmm_matrix_t;
public:
    /**
     * @brief Default constructor 
     * @details
     * Initializes all dynamically stored members to nullptr
     */
	fmm_matlab() 
		: p_surf_mesh(nullptr)
		, p_field_mesh(nullptr)
		, p_tree(nullptr)
		, p_lists(nullptr)
		, p_fmm(nullptr)
		, p_slp_matrix(nullptr)
		, p_dlp_matrix(nullptr)
		, m_accuracy(3.0)
		, m_far_field_order(5)
	{
		
	}
	
	/**
	 * @brief Create the surface mesh 
	 * @param[in] surf_nodes Surface node locations (N x 3 matrix)
	 * @param[in] surf_elems Surface element connectivity (E x 9 matrix)
	 */
	void create_surface_mesh(dmex_matrix_t const& surf_nodes, dmex_matrix_t const& surf_elems)
	{
		p_surf_mesh = new mesh_t(NiHu::create_mesh(surf_nodes, surf_elems, NiHu::quad_1_tag()));
	}
	
	/**
	 * @brief Create the field mesh
	 * @param[in] field_nodes Field node locations (N x 3 matrix)
	 * @param[in] field_elems Field element connectivity (E x 9 matrix)
	 */
	void create_field_mesh(dmex_matrix_t const& field_nodes, dmex_matrix_t const& field_elems)
    {
        p_field_mesh = new mesh_t(NiHu::create_mesh(field_nodes, field_elems, NiHu::quad_1_tag()));
    }
	
	/**
	 * @brief Create tree and interaction lists
	 * @tparam DivideDerived Cluster tree division method 
	 * @param divide Instance of the division method
	 * @details
	 * In surface mode (when no field mesh is present) the sources and receivers
	 * are the same in the cluster tree. In field mode, the sources are on the 
	 * surface mesh, while the receivers are on the field mesh.
	 */
	template <class DivideDerived>
	void create_tree(NiHu::fmm::divide_base<DivideDerived> const &divide)
	{
        if (p_field_mesh == nullptr) {
            // Create cluster tree for surface only
            p_tree = new cluster_tree_t(
                NiHu::fmm::create_elem_center_iterator(p_surf_mesh->begin<elem_t>()),
                NiHu::fmm::create_elem_center_iterator(p_surf_mesh->end<elem_t>()),
                divide
            );
        } else {
            // Create cluster tree for surface and mesh
            p_tree = new cluster_tree_t(
                NiHu::fmm::create_elem_center_iterator(p_surf_mesh->begin<elem_t>()),
                NiHu::fmm::create_elem_center_iterator(p_surf_mesh->end<elem_t>()),
                NiHu::fmm::create_elem_center_iterator(p_field_mesh->begin<elem_t>()),
                NiHu::fmm::create_elem_center_iterator(p_field_mesh->end<elem_t>()),
                divide
            );
        }
		
		// create interaction lists
		p_lists = new NiHu::fmm::interaction_lists(*p_tree);
	}
	
	/**
	 * @brief Create FMM matrices with precomputation
	 */
	void create_matrices()
    {
        if (p_field_mesh == nullptr)
            create_surface_matrices();
        else
            create_field_matrices();
    }
	
	void create_surface_matrices()
	{
		// create the FMM instance
		p_fmm = new fmm_t(m_wave_number);
		p_fmm->set_accuracy(m_accuracy);

		// initialize tree data
		p_fmm->init_level_data(*p_tree);
		for (size_t c = 0; c < p_tree->get_n_clusters(); ++c)
			(*p_tree)[c].set_p_level_data(&p_fmm->get_level_data((*p_tree)[c].get_level()));
		
		#ifdef NIHU_FMM_PARALLEL
			auto max_num_threads = omp_get_max_threads();
			for (size_t i = 0; i < p_tree->get_n_levels(); ++i)
				p_fmm->get_level_data(i).set_num_threads(max_num_threads);
		#endif
		
		// create functors
		auto int_fctr = NiHu::fmm::create_integrated_functor(
			test_field_tag_t(), trial_field_tag_t(), m_far_field_order, true);

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
		// Create the operator collection for the DLP matrix
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

		// Create the collection of additional operators for the SLP matrix
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
		auto dlp_pre_collection = dlp_cix_collection.transform(pre_fctr);
		auto slp_pre_collection = slp_cix_collection.transform(pre_fctr);
		
		// get assembly times
		m_assembly_times.fill_times(dlp_pre_collection);

		// create slp matrix object
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
		p_dlp_matrix = new fmm_matrix_t(create_fmm_matrix(
			dlp_pre_collection,
			*p_tree, *p_lists));
	}
	
	void create_field_matrices()
	{
		// create the FMM instance
		p_fmm = new fmm_t(m_wave_number);
		p_fmm->set_accuracy(m_accuracy);

		// initialize tree data
		p_fmm->init_level_data(*p_tree);
		for (size_t c = 0; c < p_tree->get_n_clusters(); ++c)
			(*p_tree)[c].set_p_level_data(&p_fmm->get_level_data((*p_tree)[c].get_level()));
		
		#ifdef NIHU_FMM_PARALLEL
			auto max_num_threads = omp_get_max_threads();
			for (size_t i = 0; i < p_tree->get_n_levels(); ++i)
				p_fmm->get_level_data(i).set_num_threads(max_num_threads);
		#endif
		
		// create functors
		auto int_fctr = NiHu::fmm::create_integrated_functor(
			test_field_tag_t(), trial_field_tag_t(), m_far_field_order, false);

		auto const &trial_space = NiHu::constant_view(*p_surf_mesh);
		auto const &test_space  = NiHu::dirac(NiHu::constant_view(*p_field_mesh));

		auto idx_fctr = NiHu::fmm::create_indexed_functor(
			test_space.template field_begin<test_field_t>(),
			test_space.template field_end<test_field_t>(),
			trial_space.template field_begin<trial_field_t>(),
			trial_space.template field_end<trial_field_t>(),
			*p_tree);

		auto pre_fctr = NiHu::fmm::create_precompute_functor(*p_tree, *p_lists);

		// Create the operator collection for the DLP matrix
		auto dlp_collection = NiHu::fmm::create_fmm_operator_collection(
			int_fctr(p_fmm->template create_p2p<0, 1>()),
			int_fctr(p_fmm->template create_p2m<1>()),
			int_fctr(p_fmm->template create_p2l<1>()),
			int_fctr(p_fmm->template create_m2p<0>()),
			int_fctr(p_fmm->template create_l2p<0>()),
			p_fmm->create_m2m(),
			p_fmm->create_m2l(),
			p_fmm->create_l2l()
		);

		// Create the collection of additional operators for the RHS matrix
		auto slp_collection = NiHu::fmm::create_fmm_operator_collection(
			int_fctr(p_fmm->template create_p2p<0, 0>()),
			int_fctr(p_fmm->template create_p2m<0>()),
			int_fctr(p_fmm->template create_p2l<0>())
		);

		// indexing
		auto dlp_cix_collection = dlp_collection.transform(idx_fctr);
		auto slp_cix_collection = slp_collection.transform(idx_fctr);

		// precomputation
		auto dlp_pre_collection = dlp_cix_collection.transform(pre_fctr);
		auto slp_pre_collection = slp_cix_collection.transform(pre_fctr);
		
		// get assembly times
		m_assembly_times.fill_times(dlp_pre_collection);
		
		// create slp matrix object
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
	void mvp_dlp(Eigen::MatrixBase<LhsDerived> &res, Eigen::MatrixBase<RhsDerived> const &src)
	{
		res = (*p_dlp_matrix) * src;
	}
	
	template <class LhsDerived, class RhsDerived>
	void mvp_slp(Eigen::MatrixBase<LhsDerived> &&res, Eigen::MatrixBase<RhsDerived> const &src)
	{
		res = (*p_slp_matrix) * src;
	}
	
	template <class LhsDerived, class RhsDerived>
	void mvp_slp(Eigen::MatrixBase<LhsDerived> &res, Eigen::MatrixBase<RhsDerived> const &src)
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
	
	void set_far_field_order(size_t far_field_order)
    {
        m_far_field_order = far_field_order;
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
	
	size_t get_rows() const
	{
		return p_slp_matrix->rows();
	}
	
	size_t get_cols() const
	{
		return p_slp_matrix->cols();
	}
	
	mesh_t const * get_surface_mesh() const
	{
		return p_surf_mesh;
	}
	
	mesh_t const * get_field_mesh() const
	{
		return p_field_mesh;
	}
	
	fmm_assembly_times const &get_assembly_times() const
	{
		return m_assembly_times;
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
    size_t m_far_field_order;
	
	fmm_assembly_times m_assembly_times;
};

fmm_matlab *p = nullptr;

void usage(int nrhs, mxArray const *prhs[])
{
	if (nrhs < 2) {
		// Print general usage 
		mexPrintf(NIHU_THIS_MEX_NAME" -- general usage:\n\n");
		mexPrintf("Use this MEX function as " NIHU_THIS_MEX_NAME "(command), \n"
			"where the command string can be the following.\n\n");
		mexPrintf("  'help'         Display help information.\n");
		mexPrintf("  'cleanup'      Clean up allocated structures.\n");
		mexPrintf("  'init'         Initialize the FMM object for further operations.\n");
		mexPrintf("  'matrix'       Assemble and FMM matrices with precomputing.\n");
		mexPrintf("  'mesh'         Initialize surface and/or field meshes.\n");
		mexPrintf("  'mvp_dlp'      Compute matrix vector product using the DLP matrix.\n");
		mexPrintf("  'mvp_slp'      Compute matrix vector product using the SLP matrix.\n");
		mexPrintf("  'set'          Set values for parameters.\n");
		mexPrintf("  'print_times'  Print matrix assembly times.\n");
		mexPrintf("  'print_tree'   Print tree information.\n");
		mexPrintf("  'tree'         Build the cluster tree.\n");
		mexPrintf("\n");
		mexPrintf("For further help on specific commands use:\n");
		mexPrintf("  " NIHU_THIS_MEX_NAME "('help', command)\n");
	} else {
		if (mxIsChar(prhs[1]) != 1) {
			mexErrMsgIdAndTxt("NiHu:" NIHU_THIS_MEX_NAME ":invalid_input",
				"Parameter name must be a string for the command \"help\".");
		}
		char const *help_cmd = mxArrayToString(prhs[1]);
		
		// Help for 'init'
		if (!strcmp(help_cmd, "init")) {
			mexPrintf(NIHU_THIS_MEX_NAME "('init')\n\n");
			mexPrintf("Initializes the data structures needed for the FMBEM computation.\n\n");
			mexPrintf("This command does not take additional arguments. Successive\n"
				      "calls to the command 'init' are refused.\n\n");
			mexPrintf("See also the command 'cleanup'.\n\n");
		}
		
		// Help for 'cleanup'
		else if (!strcmp(help_cmd, "cleanup")) {
			mexPrintf(NIHU_THIS_MEX_NAME "('cleanup')\n\n");
			mexPrintf("Cleans up the allocated structures of the FMBEM computation.\n\n");
			mexPrintf("This command does not take additional arguments.\n"
			          "The 'cleanup' command can only be called after the 'init'\n"
					  "command was successfully performed.\n"
			          "Successive calls to the command 'cleanup' are refused.\n\n");
			mexPrintf("See also the command 'init'.\n\n");
		}

		// Help for 'set'
		else if (!strcmp(help_cmd, "set")) {
			mexPrintf(NIHU_THIS_MEX_NAME "('set', 'param1', value1, 'param2', value2, ...)\n\n");
			mexPrintf("Sets the values to the given parameters.\n\n");
			mexPrintf("The parameters 'param1', 'param2', ... are the names of the parameters\n"
			          "to be set. Each name is followed by the corresponding value.\n"
					  "The following parameters can be set:\n\n");
			mexPrintf("  'accuracy'         Scalar parameter C\n"
					  "        Controls the accuracy of the series expansion of the Green's function.\n"
					  "        The resulting expansion length L is found from the wave number k, the\n"
					  "        cluster diameter d and the accuracy parameter C as:\n"
					  "        L = ceil(kd + C * log(kd + pi))\n"
					  "        Default value: C = 3.\n\n");
			mexPrintf("  'far_field_order'  Positive integer far field quadrature order\n"
					  "        Controls the number of points of the Gaussian quadrature used in far\n"
					  "        field integration. Default value: far_field_order = 5.\n\n");
			mexPrintf("  'wave_number'  Scalar wave number k\n"
					  "        The wave number parmeter k for the Helmholtz problem.\n\n");
		}
		
		// Help for 'mesh'
		else if (!strcmp(help_cmd, "mesh")) {
			mexPrintf(NIHU_THIS_MEX_NAME "('mesh', s_nodes, s_elems)\n\n");
			mexPrintf("Sets up the surface mesh for the FMBEM computation.\n");
			mexPrintf("The matrices s_nodes and s_elems contain the nodal locations and\n"
			          "the element connectivity of the surface mesh. This function supports\n"
					  "quadrilateral elements only.\n");
			mexPrintf("In surface mode the SLP and DLP matrices are evaluted using the\n"
			          "Burton-Miller formulation for avoiding fictitious eigenfrequencies.\n"
					  "(Nearly-)Singular integrals are treated using special techniques.\n"
					  "The SLP and DLP matrices form the following relation between the\n"
					  "surface pressure ps and its normal derivative qs:\n"
					  "DLP * ps = SLP * qs\n\n");
			mexPrintf(NIHU_THIS_MEX_NAME "('mesh', s_nodes, s_elems, f_nodes, f_elems)\n\n");
			mexPrintf("Sets up the surface and field meshes for the FMBEM computation.\n");
			mexPrintf("The matrices s_nodes, s_elems, f_nodes, f_elems contain the nodal\n" 
					  "locations and the element connectivity of the surface and field meshes,\n"
					  "respectively. This function supports quadrilateral elements only.\n");
			mexPrintf("In field mode, no singular treatment is performed, and the\n"
					  "simple collocational form is used. The field point pressure pf\n"
					  "is found then by the matrix-vector products:\n"
					  "pf = DLP * ps - SLP * qs\n\n");
			mexPrintf("The 'mesh' command must be called after the 'init' command was\n"
			          "successfully completed.\n\n");
			mexPrintf("See also the commands 'matrix, 'mvp_dlp', and 'mvp_slp'.\n");
			mexPrintf("See also the Matlab function extract_core_mesh\n\n");
		}
		
		// Help for 'tree'
		else if (!strcmp(help_cmd, "tree")) {
			mexPrintf(NIHU_THIS_MEX_NAME "('tree', option, param)\n\n");
			mexPrintf("Creates the cluster tree and assembles the interaction lists.\n");
			mexPrintf("The option parameter controls the tree building strategy, with"
			          "the following possible configurations:\n\n");
			mexPrintf("  'divide_depth'      Parameter: positive integer L.\n"
					  "        Divides the cluster tree such that the leaf level is level L.\n\n");
			mexPrintf("  'divide_diameter'   Parameter: positive real D.\n"
			          "        Divides each cluster as long as its diameter is greater than D.\n\n");
			mexPrintf("  'divide_num_nodes'  Parameter: position integer N.\n"
					  "        Divides each cluster as long as it has more nodes than N.\n\n");
			mexPrintf("See also the command 'print_tree'\n\n");
		}
		
		// Help for 'matrix'
		else if (!strcmp(help_cmd, "matrix")) {
			mexPrintf(NIHU_THIS_MEX_NAME "('matrix')\n\n");
			mexPrintf("Assembles the FMBEM matrices for fast matrix-vector products.\n"
			          "The SLP and DLP matrices can then be used for solving the the\n"
					  "BEM system in mesh mode, or computing the field point pressure\n"
					  "in field mode.\n\n");
			mexPrintf("See also the commands 'mesh', 'mvp_dlp', and 'mvp_slp'.\n\n");
		}
		
		else if(!strcmp(help_cmd, "mvp_slp") || !strcmp(help_cmd, "mvp_dlp")) {
			mexPrintf("r = " NIHU_THIS_MEX_NAME "('mvp_dlp', ps)\n");
			mexPrintf("r = " NIHU_THIS_MEX_NAME "('mvp_slp', qs)\n\n");
			mexPrintf("Evaluate fast matrix-vector product using the SLP or DLP matrix.\n");
			mexPrintf("The input vector ps or qs must have the same number of elements\n"
			          "as the number of elements in the surface mesh. In surface mode,\n"
					  "the result vector r will have the same number of elements as the\n"
					  "surface mesh, while in field mode the vector r will have the same\n"
					  "number of elements as the field mesh.\n\n");
		}
		
		// Diagnostic commands and help
		else if(!strcmp(help_cmd, "help") || !strcmp(help_cmd, "print_tree") || !strcmp(help_cmd, "print_times")) {
			mexPrintf(NIHU_THIS_MEX_NAME "('help', '%s')\n\n", help_cmd);
			mexPrintf("There is no further documentation for this command.\n\n");
		}
		
		// Unsupported option
		else {
			mexErrMsgIdAndTxt("NiHu:" NIHU_THIS_MEX_NAME ":unknown_option",
				"The option '%s' is not recognized for the command \"help\".", help_cmd);
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
	char const *input_option;
	
	/* input must be a string */
    if ( mxIsChar(prhs[0]) != 1) {
		mexErrMsgIdAndTxt("NiHu:" NIHU_THIS_MEX_NAME ":invalid_input",
              "First input parameter must be a string, use " NIHU_THIS_MEX_NAME "('help') to see usage");
	}
	input_option = mxArrayToString(prhs[0]);
	
	if (!strcmp(input_option, "help")) {
		usage(nrhs, prhs);
	}
	
	// Init option
	else if (!strcmp(input_option, "init")) {
		if (p != nullptr) {
			mexErrMsgIdAndTxt("NiHu:" NIHU_THIS_MEX_NAME ":invalid_state",
				"Command \"%s\" called in invalid state", input_option);
		} else {
			p = new fmm_matlab;
		}
	}
	
	// Set option
	else if (!strcmp(input_option, "set")) {
		if (p == nullptr) {
			mexErrMsgIdAndTxt("NiHu:" NIHU_THIS_MEX_NAME ":invalid_state",
				"Command \"%s\" called in invalid state", input_option);
		}
		
		// Go through parameter pairs
		for (int i = 0; i < (nrhs - 1) / 2; ++i)
		{
			if ( mxIsChar(prhs[2*i + 1]) != 1) {
				mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_surf:invalid_input",
				"Parameter name must be a string for the command \"set\".");
			}
			
			char const *what_to_set = mxArrayToString(prhs[2 * i + 1]);
			if (!strcmp(what_to_set, "wave_number")) {
				double k = mxGetScalar(prhs[2 * i + 2]);
				p->set_wave_number(k);
			} else if (!strcmp(what_to_set, "accuracy")) {
				double accuracy = mxGetScalar(prhs[2 * i + 2]);
				p->set_accuracy(accuracy);
			} else if (!strcmp(what_to_set, "far_field_order")) {
                size_t far_field_order = size_t(mxGetScalar(prhs[2 * i + 2]));
                p->set_far_field_order(far_field_order);
            } else {
				mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_surf:invalid_parameter",
					"Unknown parameter name to set: \"%s\"", what_to_set);
			}
		}
	}
	
	// Mesh option
	else if (!strcmp(input_option, "mesh"))
	{
		if (p == nullptr) {
			mexErrMsgIdAndTxt("NiHu: " NIHU_THIS_MEX_NAME ":invalid_state",
				"Command \"%s\" called in invalid state", input_option);
		} else {
            if (nrhs == 3) {
                p->create_surface_mesh(dmex_matrix_t(prhs[1]), dmex_matrix_t(prhs[2]));
            } else if (nrhs == 5) {
                p->create_surface_mesh(dmex_matrix_t(prhs[1]), dmex_matrix_t(prhs[2]));
                p->create_field_mesh(dmex_matrix_t(prhs[3]), dmex_matrix_t(prhs[4]));
            } else {
                mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_mex:invalid_input",
                    "Command \"mesh\" called with invalid input");
            }
		}
	}
	
	else if (!strcmp(input_option, "tree"))
	{
		if (p == nullptr || p->get_surface_mesh() == nullptr) {
			mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_mex:invalid_state",
				"Command \"tree\" called in invalid state");
		} else {
			if ( mxIsChar(prhs[1]) != 1) {
				mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_surf:invalid_input",
				"Division method name must be a string for the command \"%s\".", input_option);
			}
			char const *divide_option = mxArrayToString(prhs[1]);
			if (!strcmp(divide_option, "divide_depth")) {
				p->create_tree(NiHu::fmm::divide_depth(size_t(mxGetScalar(prhs[2]))));
			} else if (!strcmp(divide_option, "divide_num_nodes")) {
				p->create_tree(NiHu::fmm::divide_num_nodes(size_t(mxGetScalar(prhs[2]))));
			} else if (!strcmp(divide_option, "divide_diameter")) {
				p->create_tree(NiHu::fmm::divide_diameter(mxGetScalar(prhs[2])));
			} else {
				mexErrMsgIdAndTxt("NiHu:" NIHU_THIS_MEX_NAME ":invalid_divide_option",
					"Unknown divide option: \"%s\"", divide_option);
			}
		}
	}
	
	// Matrix assembly
	else if (!strcmp(input_option, "matrix"))
	{
		if (p == nullptr)
		{
			mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_mex:invalid_state",
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
			mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_mex:invalid_state",
				"Command \"%s\" called in invalid state", input_option);
		}
		else
		{
			cmex_matrix_t xct(prhs[1]);
			cmex_matrix_t res(p->get_rows(), 1, plhs[0]);
			#if MX_HAS_INTERLEAVED_COMPLEX
				p->mvp_dlp(res.col(0), xct.col(0));
			#else
				// Copying is needed in case of old Matlab complex storage
				cmatrix_t x(xct.rows(), xct.cols());
				cmatrix_t r(p->get_rows(), res.cols());
				for (int i = 0; i < xct.rows(); ++i)
					x(i, 0) = xct(i, 0);
				
				p->mvp_dlp(r.col(0), x.col(0));
				
				for (int i = 0; i < r.rows(); ++i)
					res(i, 0) = r(i, 0);
			#endif
		}
	}
	
	else if (!strcmp(input_option, "mvp_slp"))
	{
		if (p == nullptr)
		{
			mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_mex:invalid_state",
				"Command \"%s\" called in invalid state", input_option);
		}
		else
		{
			cmex_matrix_t xct(prhs[1]);
			cmex_matrix_t res(p->get_rows(), 1, plhs[0]);
			#if MX_HAS_INTERLEAVED_COMPLEX
				p->mvp_slp(res.col(0), xct.col(0));
			#else
				// Copying is needed in case of old Matlab complex storage
				cmatrix_t x(xct.rows(), xct.cols());
				cmatrix_t r(p->get_rows(), res.cols());
				for (int i = 0; i < xct.rows(); ++i)
					x(i, 0) = xct(i, 0);
				
				p->mvp_slp(r.col(0), x.col(0));
				
				for (int i = 0; i < r.rows(); ++i)
					res(i, 0) = r(i, 0);
				
			#endif
		}
	}
	
	// Cleanup option 
	else if (!strcmp(input_option, "cleanup"))
	{
		if (p == nullptr) {
			mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_mex:invalid_state",
				"Command \"%s\" called in invalid state", input_option);
		} else {
			delete p;
			p = nullptr;
		}
	}
	
	// Diagnostics
	// Print matrix assembly times
	else if (!strcmp(input_option, "print_times")) {
		if (p == nullptr) {
			mexErrMsgIdAndTxt("NiHu:" NIHU_THIS_MEX_NAME ":invalid_state",
				"Command \"%s\" called in invalid state", input_option);
		} else {
			fmm_assembly_times const& times = p->get_assembly_times();
			mexPrintf("Assembly times (in microseconds):\n");
			mexPrintf("\tM2M: %10llu\n", times.get_time(NiHu::fmm::fmm_timer::M2M));
			mexPrintf("\tM2L: %10llu\n", times.get_time(NiHu::fmm::fmm_timer::M2L));
			mexPrintf("\tL2L: %10llu\n", times.get_time(NiHu::fmm::fmm_timer::L2L));
			mexPrintf("\tP2M: %10llu\n", times.get_time(NiHu::fmm::fmm_timer::P2M));
			mexPrintf("\tP2L: %10llu\n", times.get_time(NiHu::fmm::fmm_timer::P2L));
			mexPrintf("\tM2P: %10llu\n", times.get_time(NiHu::fmm::fmm_timer::M2P));
			mexPrintf("\tL2P: %10llu\n", times.get_time(NiHu::fmm::fmm_timer::L2P));
			mexPrintf("\tP2P: %10llu\n", times.get_time(NiHu::fmm::fmm_timer::P2P));
		}
	}
	
	// Print tree structure
	else if (!strcmp(input_option, "print_tree")) {
		if (p == nullptr) {
			mexErrMsgIdAndTxt("NiHu:" NIHU_THIS_MEX_NAME ":invalid_state",
				"Command \"%s\" called in invalid state", input_option);
		} else {
			p->print_tree();
		}
		
	}
	
	// The option was not valid
	else {
		mexErrMsgIdAndTxt("NiHu:helmholtz_3d_hf_fmm_mex:invalid_option",
			"Unknown input option: \"%s\"", input_option);
	}	

}
