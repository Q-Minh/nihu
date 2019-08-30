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

typedef NiHu::covariance_kernel<NiHu::space_2d<> > kernel_t;
typedef NiHu::fmm::black_box_fmm<kernel_t> bbfmm_t;
typedef bbfmm_t::cluster_t cluster_t;
typedef NiHu::fmm::cluster_tree<cluster_t> cluster_tree_t;
typedef NiHu::quad_1_volume_elem elem_t;
typedef NiHu::field_view<elem_t, NiHu::field_option::constant> field_t;
typedef NiHu::type2tag<field_t>::type field_tag_t;

int main(int argc, char *argv[])
{
	double sigma = 1.;
	double length = 1.;
	kernel_t kernel(sigma, length);

	bbfmm_t bbfmm(kernel);

	std::string meshname(argv[1]);
	auto mesh = NiHu::read_off_mesh(meshname, NiHu::quad_1_volume_tag());
	auto space = NiHu::constant_view(mesh);

	cluster_tree_t tree(
		NiHu::fmm::create_elem_center_iterator(mesh.template begin<elem_t>()),
		NiHu::fmm::create_elem_center_iterator(mesh.template end<elem_t>()),
		NiHu::fmm::divide_num_nodes(10)
	);
	for (size_t i = 0; i < tree.get_n_clusters(); ++i)
		tree[i].set_chebyshev_order(5);

	NiHu::fmm::interaction_lists lists(tree);

	size_t quadrature_order = 5;
	auto int_fctr = NiHu::fmm::create_integrated_functor(
		field_tag_t(), field_tag_t(), quadrature_order, true);

	auto idx_fctr = NiHu::fmm::create_indexed_functor(
		space.field_begin<field_t>(),
		space.field_end<field_t>(),
		space.field_begin<field_t>(),
		space.field_end<field_t>(),
		tree
	);

	auto pre_fctr = NiHu::fmm::create_precompute_functor(tree, lists);

	auto collection = NiHu::fmm::create_fmm_operator_collection(
		int_fctr(bbfmm.create_p2p()),
		int_fctr(bbfmm.create_p2m()),
		int_fctr(bbfmm.create_p2l()),
		bbfmm.create_m2m(),
		bbfmm.create_m2l(),
		bbfmm.create_l2l(),
		int_fctr(bbfmm.create_l2p()),
		int_fctr(bbfmm.create_m2p())
	);

	auto pre_collection = collection.transform(idx_fctr).transform(pre_fctr);

	auto mat = NiHu::fmm::create_fmm_matrix(pre_collection, tree, lists);

	Eigen::Matrix<double, Eigen::Dynamic, 1> xct(mat.cols(), 1);
	xct.setConstant(1.0);
	auto res = mat * xct;

	std::cout << res << std::endl;

	return 0;
}

