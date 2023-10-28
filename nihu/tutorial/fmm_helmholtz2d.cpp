#include "nihu/fmm/helmholtz_2d_wb_fmm.hpp"
#include "nihu/fmm/lists.hpp"
#include "nihu/fmm/fmm_indexed.hpp"
#include "nihu/fmm/fmm_precompute.hpp"
#include "nihu/fmm/fmm_matrix.hpp"

//! [typedefs]
typedef double wave_number_t;
typedef NiHu::fmm::helmholtz_2d_wb_fmm<wave_number_t> fmm_t;
typedef fmm_t::cluster_t cluster_t;
typedef NiHu::fmm::cluster_tree<cluster_t> cluster_tree_t;
typedef cluster_t::location_t location_t;
//! [typedefs]

int main()
{
	//! [source and receiver nodes]
	size_t N_src = 5000;	// number of source nodes
	location_t *src_begin = new location_t[N_src];
	for (size_t i = 0; i < N_src; ++i)
		src_begin[i].setRandom();
	size_t N_rec = 3000;	// number of receiver nodes
	location_t *rec_begin = new location_t[N_rec];
	for (size_t i = 0; i < N_rec; ++i)
	{
		rec_begin[i].setRandom();
		rec_begin[i](0) += 2.0;
	}
	//! [source and receiver nodes]

	//! [tree and lists]
	cluster_tree_t tree(src_begin, src_begin + N_src,
		rec_begin, rec_begin + N_rec,
		NiHu::fmm::divide_num_nodes(10));
	
	NiHu::fmm::interaction_lists lists(tree);
	//! [tree and lists]

	//! [fmm initialization]
	wave_number_t k = 2.0;
	fmm_t fmm(k);

	fmm.init_level_data(tree);

	for (size_t c = 0; c < tree.get_n_clusters(); ++c)
		tree[c].set_p_level_data(&fmm.get_level_data(tree[c].get_level()));
	//! [fmm initialization]

	//! [operator manipulations]
	auto ops = NiHu::fmm::create_fmm_operator_collection(
		fmm.create_p2p<0, 0>(),
		fmm.create_p2m<0>(),
		fmm.create_p2l<0>(),
		fmm.create_m2p<0>(),
		fmm.create_l2p<0>(),
		fmm.create_m2m(),
		fmm.create_l2l(),
		fmm.create_m2l()
	);

	auto p2p = fmm.create_p2p<0, 0>();
	std::cout << p2p(*rec_begin, *src_begin) << std::endl;

	auto idx_fctr = NiHu::fmm::create_indexed_functor(
		src_begin, src_begin + N_src,
		rec_begin, rec_begin + N_rec,
		tree);

	auto pre_fctr = NiHu::fmm::create_precompute_functor(tree, lists);

	auto pre_collection = ops.transform(idx_fctr).transform(pre_fctr);
	//! [operator manipulations]

	//! [matrix]
	auto mat = NiHu::fmm::create_fmm_matrix(pre_collection, tree, lists);

	Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> xct(N_src, 1);
	xct.setConstant(1.0);

	auto resp = mat * xct;
	//! [matrix]

	delete[] src_begin;
	delete[] rec_begin;

	return 0;
}
