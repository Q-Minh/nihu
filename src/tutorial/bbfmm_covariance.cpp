#include "library/covariance_kernel.hpp"
#include "fmm/black_box_fmm.hpp"
#include "fmm/cluster_tree.hpp"
#include "fmm/elem_center_iterator.hpp"
#include "fmm/fmm_operator_collection.hpp"
#include "fmm/fmm_integrated.hpp"
#include "fmm/fmm_indexed.hpp"
#include "fmm/fmm_precompute.hpp"
#include "fmm/fmm_matrix.hpp"
#include "fmm/lists.hpp"
#include "interface/read_off_mesh.hpp"
#include "library/lib_element.hpp"
#include "core/mesh.hpp"
#include "core/function_space.hpp"
#include "core/field.hpp"

typedef NiHu::covariance_kernel<NiHu::space_2d<> > kernel_t;
typedef NiHu::quad_1_volume_elem elem_t;
typedef NiHu::type2tag<elem_t>::type elem_tag_t;
typedef NiHu::field_view<elem_t, NiHu::field_option::constant> field_t;
typedef NiHu::type2tag<field_t>::type field_tag_t;
typedef NiHu::fmm::black_box_fmm<kernel_t> fmm_t;
typedef fmm_t::cluster_t cluster_t;
typedef NiHu::fmm::cluster_tree<cluster_t> cluster_tree_t;

typedef NiHu::mesh<tmp::vector<elem_t> > mesh_t;
typedef NiHu::function_space_view<mesh_t, NiHu::field_option::constant> space_t;

class C
{

public:
	C(int argc, char const *argv[])
		: m_kernel(1., .1)
		, m_bbfmm(m_kernel)
		, m_mesh(NiHu::read_off_mesh(argv[1], elem_tag_t()))
		, m_space(NiHu::constant_view(m_mesh))
		, m_tree(NiHu::fmm::create_elem_center_iterator(m_mesh.begin<elem_t>()),
			NiHu::fmm::create_elem_center_iterator(m_mesh.end<elem_t>()),
			NiHu::fmm::divide_num_nodes(10))
		, m_lists(m_tree)
		, m_int_fctr(NiHu::fmm::create_integrated_functor(
			field_tag_t(), field_tag_t(), 5, true))
		, m_idx_fctr(NiHu::fmm::create_indexed_functor(
			m_space.field_begin<field_t>(),
			m_space.field_end<field_t>(),
			m_space.field_begin<field_t>(),
			m_space.field_end<field_t>(),
			m_tree))
		, m_pre_fctr(NiHu::fmm::create_precompute_functor(m_tree, m_lists))
	{
		for (size_t i = 0; i < m_tree.get_n_clusters(); ++i)
			m_tree[i].set_chebyshev_order(5);

		m_mat = new mat_t(NiHu::fmm::create_fmm_matrix(
			m_pre_fctr(m_idx_fctr(m_int_fctr(m_bbfmm.create_p2p()))),
			m_pre_fctr(m_idx_fctr(m_int_fctr(m_bbfmm.create_p2m()))),
			m_pre_fctr(m_idx_fctr(m_int_fctr(m_bbfmm.create_p2l()))),
			m_pre_fctr(m_idx_fctr(m_int_fctr(m_bbfmm.create_m2p()))),
			m_pre_fctr(m_idx_fctr(m_int_fctr(m_bbfmm.create_l2p()))),
			m_pre_fctr(m_idx_fctr(m_bbfmm.create_m2m())),
			m_pre_fctr(m_idx_fctr(m_bbfmm.create_l2l())),
			m_pre_fctr(m_idx_fctr(m_bbfmm.create_m2l())),
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
	decltype(NiHu::fmm::create_integrated_functor(
		field_tag_t(), field_tag_t(), 5, true)) const m_int_fctr;
	decltype(NiHu::fmm::create_indexed_functor(
		m_space.field_begin<field_t>(),
		m_space.field_end<field_t>(),
		m_space.field_begin<field_t>(),
		m_space.field_end<field_t>(),
		m_tree)) const m_idx_fctr;
	decltype(NiHu::fmm::create_precompute_functor(m_tree, m_lists)) const m_pre_fctr;

	typedef decltype(NiHu::fmm::create_fmm_matrix(
		m_pre_fctr(m_idx_fctr(m_int_fctr(m_bbfmm.create_p2p()))),
		m_pre_fctr(m_idx_fctr(m_int_fctr(m_bbfmm.create_p2m()))),
		m_pre_fctr(m_idx_fctr(m_int_fctr(m_bbfmm.create_p2l()))),
		m_pre_fctr(m_idx_fctr(m_int_fctr(m_bbfmm.create_m2p()))),
		m_pre_fctr(m_idx_fctr(m_int_fctr(m_bbfmm.create_l2p()))),
		m_pre_fctr(m_idx_fctr(m_bbfmm.create_m2m())),
		m_pre_fctr(m_idx_fctr(m_bbfmm.create_l2l())),
		m_pre_fctr(m_idx_fctr(m_bbfmm.create_m2l())),
		m_tree,
		m_lists)) mat_t;

	mat_t *m_mat;
};

int main(int argc, char const *argv[])
{
	C c(argc, argv);
	return 0;
}

