#include "nihu/library/laplace_kernel.hpp"
#include "nihu/library/laplace_singular_integrals.hpp"
#include "nihu/library/laplace_nearly_singular_integrals.hpp"
#include "nihu/fmm/black_box_fmm.hpp"
#include "nihu/fmm/cluster_tree.hpp"
#include "nihu/fmm/elem_center_iterator.hpp"
#include "nihu/fmm/fmm_operator_collection.hpp"
#include "nihu/fmm/fmm_integrated.hpp"
#include "nihu/fmm/fmm_indexed.hpp"
#include "nihu/fmm/fmm_precompute.hpp"
#include "nihu/fmm/fmm_matrix.hpp"
#include "nihu/fmm/lists.hpp"
#include "interface/read_off_mesh.hpp"
#include "nihu/library/lib_element.hpp"
#include "nihu/core/mesh.hpp"
#include "nihu/core/function_space.hpp"
#include "nihu/core/field.hpp"

typedef NiHu::laplace_3d_DLP_kernel kernel_t;
typedef NiHu::quad_1_elem elem_t;
typedef NiHu::type2tag<elem_t>::type elem_tag_t;
typedef NiHu::field_view<elem_t, NiHu::field_option::constant> field_t;
typedef NiHu::type2tag<field_t>::type field_tag_t;
typedef NiHu::fmm::black_box_fmm<kernel_t> fmm_t;
typedef fmm_t::cluster_t cluster_t;
typedef NiHu::fmm::cluster_tree<cluster_t> cluster_tree_t;

typedef NiHu::mesh<tmp::vector<elem_t> > mesh_t;
typedef NiHu::function_space_view<mesh_t, NiHu::field_option::constant> space_t;

int main(int argc, char const *argv[])
{
	kernel_t kernel;
	fmm_t bbfmm(kernel);
	auto mesh = NiHu::read_off_mesh(argv[1], elem_tag_t());
	auto const &space = NiHu::constant_view(mesh);
	cluster_tree_t tree(NiHu::fmm::create_elem_center_iterator(mesh.begin<elem_t>()),
		NiHu::fmm::create_elem_center_iterator(mesh.end<elem_t>()),
		NiHu::fmm::divide_num_nodes(10));
	NiHu::fmm::interaction_lists lists(tree);
	auto int_fctr = NiHu::fmm::create_integrated_functor(
		field_tag_t(), field_tag_t(), 5, true);
	auto idx_fctr = NiHu::fmm::create_indexed_functor(
		space.field_begin<field_t>(),
		space.field_end<field_t>(),
		space.field_begin<field_t>(),
		space.field_end<field_t>(),
		tree);
	auto pre_fctr = NiHu::fmm::create_precompute_functor(tree, lists);

	for (size_t i = 0; i < tree.get_n_clusters(); ++i)
		tree[i].set_chebyshev_order(5);

	auto mat = NiHu::fmm::create_fmm_matrix(
		pre_fctr(idx_fctr(int_fctr(bbfmm.create_p2p()))),
		pre_fctr(idx_fctr(int_fctr(bbfmm.create_p2m()))),
		pre_fctr(idx_fctr(int_fctr(bbfmm.create_p2l()))),
		pre_fctr(idx_fctr(int_fctr(bbfmm.create_m2p()))),
		pre_fctr(idx_fctr(int_fctr(bbfmm.create_l2p()))),
		pre_fctr(idx_fctr(bbfmm.create_m2m())),
		pre_fctr(idx_fctr(bbfmm.create_l2l())),
		pre_fctr(idx_fctr(bbfmm.create_m2l())),
		tree,
		lists);

	Eigen::Matrix<double, Eigen::Dynamic, 1> xct(mat.cols(), 1);
	xct.setConstant(1.0);
	auto res = mat * xct;

	std::cout << res << std::endl;
}


