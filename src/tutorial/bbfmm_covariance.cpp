#include "fmm/black_box_fmm.hpp"
#include "fmm/cluster_tree.hpp"
#include "fmm/divide.hpp"
#include "fmm/elem_center_iterator.hpp"
#include "fmm/lists.hpp"
#include "fmm/fmm_integrated.hpp"
#include "fmm/fmm_indexed.hpp"
#include "fmm/fmm_precompute.hpp"
#include "fmm/fmm_operator_collection.hpp"
#include "fmm/fmm_matrix.hpp"
#include "interface/read_off_mesh.hpp"
#include "library/covariance_kernel.hpp"
#include "library/lib_element.hpp"
#include "core/function_space.hpp"

typedef NiHu::exponential_covariance_kernel<NiHu::space_2d<>, NiHu::field_dimension::_1d> kernel_t;
typedef NiHu::fmm::black_box_fmm<kernel_t> fmm_t;
typedef fmm_t::cluster_t cluster_t;
typedef NiHu::fmm::cluster_tree<cluster_t> cluster_tree_t;
typedef NiHu::quad_1_volume_elem elem_t;
typedef NiHu::type2tag<elem_t> elem_tag_t;
typedef NiHu::mesh<tmp::vector<elem_t> > mesh_t;
typedef NiHu::function_space_view<mesh_t, NiHu::field_option::constant> space_t;
typedef NiHu::field_view<elem_t, NiHu::field_option::constant> field_t;
typedef NiHu::type2tag<field_t>::type field_tag_t;

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
	m2l_pre_t> mat_t;

class C
{
public:
	C(int argc, char const *argv[])
		: m_kernel(kernel_t::field_variance_t(1.0), .1)
		, m_bbfmm(m_kernel)
		, m_mesh(NiHu::read_off_mesh(argv[1], elem_tag_t()))
		, m_space(NiHu::constant_view(m_mesh))
		, m_tree(NiHu::fmm::create_elem_center_iterator(m_mesh.begin<elem_t>()),
			NiHu::fmm::create_elem_center_iterator(m_mesh.end<elem_t>()),
			NiHu::fmm::divide_num_nodes(10))
		, m_lists(m_tree)
	{
		for (size_t i = 0; i < m_tree.get_n_clusters(); ++i)
			m_tree[i].set_chebyshev_order(5);

		auto int_fctr = NiHu::fmm::create_integrated_functor(
			field_tag_t(), field_tag_t(), 5, true);

		auto idx_fctr = NiHu::fmm::create_indexed_functor(
			m_space.field_begin<field_t>(),
			m_space.field_end<field_t>(),
			m_space.field_begin<field_t>(),
			m_space.field_end<field_t>(),
			m_tree);

		auto pre_fctr = NiHu::fmm::create_precompute_functor(m_tree, m_lists);

		m_mat = new mat_t(NiHu::fmm::create_fmm_matrix(
			pre_fctr(idx_fctr(int_fctr(m_bbfmm.create_p2p()))),
			pre_fctr(idx_fctr(int_fctr(m_bbfmm.create_p2m()))),
			pre_fctr(idx_fctr(int_fctr(m_bbfmm.create_p2l()))),
			pre_fctr(idx_fctr(int_fctr(m_bbfmm.create_m2p()))),
			pre_fctr(idx_fctr(int_fctr(m_bbfmm.create_l2p()))),
			pre_fctr(idx_fctr(m_bbfmm.create_m2m())),
			pre_fctr(idx_fctr(m_bbfmm.create_l2l())),
			pre_fctr(idx_fctr(m_bbfmm.create_m2l())),
			m_tree,
			m_lists));

		Eigen::Matrix<double, Eigen::Dynamic, 1> xct(m_mat->cols(), 1);
		xct.setConstant(1.0);
		auto res = (*m_mat) * xct;

		delete m_mat;

		std::cout << res << std::endl;
	}

private:

	kernel_t const m_kernel;
	fmm_t const m_bbfmm;
	mesh_t const m_mesh;
	space_t const &m_space;
	cluster_tree_t m_tree;
	NiHu::fmm::interaction_lists const m_lists;
	mat_t *m_mat;
};

int main(int argc, char const *argv[])
{
	C c(argc, argv);
	return 0;
}

