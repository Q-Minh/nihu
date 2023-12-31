#include "nihu/fmm/cluster_tree.hpp"
#include "nihu/fmm/elem_center_iterator.hpp"
#include "nihu/fmm/empty_cluster.hpp"
#include "nihu/fmm/divide.hpp"
#include "interface/read_off_mesh.hpp"
#include "nihu/library/lib_element.hpp"

#include <iostream>

//! [cluster type]
static size_t const dimension = 3;
typedef NiHu::fmm::empty_cluster<dimension> cluster_t;
typedef NiHu::fmm::cluster_tree<cluster_t> cluster_tree_t;
//! [cluster type]

//! [dfs traverse]
void dfs(cluster_tree_t const &tree, size_t idx)
{
	for (auto c : tree[idx].get_children())
		dfs(tree, c);
}
//! [dfs traverse]



int main()
{
	//! [nodal locations]
	typedef cluster_t::location_t location_t;
	size_t N = 5000;	// number of nodes
	location_t *begin = new location_t[N];
	for (size_t i = 0; i < N; ++i)
		begin[i].setRandom();
	//! [nodal locations]

	//! [construct tree]
	cluster_tree_t tree(begin, begin + N, NiHu::fmm::divide_num_nodes(10));
	//! [construct tree]

	//! [construct tree alternatives]
	cluster_tree_t tree_depth(begin, begin + N, NiHu::fmm::divide_depth(6));
	cluster_tree_t tree_diam(begin, begin + N, NiHu::fmm::divide_diameter(1e-1));
	//! [construct tree alternatives]


	//! [construct tree source receiver]
	size_t N_src = 1000;
	location_t *src_begin = new location_t[N_src];
	for (size_t i = 0; i < N_src; ++i)
		src_begin[i].setRandom();

	size_t N_rec = 200;
	location_t *rec_begin = new location_t[N_rec];
	for (size_t i = 0; i < N_rec; ++i)
		rec_begin[i].setRandom();

	cluster_tree_t src_rec_tree(src_begin, src_begin + N_src,
		rec_begin, rec_begin + N_rec,
		NiHu::fmm::divide_num_nodes(10));
	//! [construct tree source receiver]

	//! [construct tree mesh]
	auto mesh = NiHu::read_off_mesh("mesh.off", NiHu::tria_1_tag());
	cluster_tree_t tree_mesh(
		NiHu::fmm::create_elem_center_iterator(mesh.template begin<NiHu::tria_1_elem>()),
		NiHu::fmm::create_elem_center_iterator(mesh.template end<NiHu::tria_1_elem>()),
		NiHu::fmm::divide_num_nodes(10));
	//! [construct tree mesh]

	//! [vector traversing]
	for (size_t i = 0; i < tree.get_n_clusters(); ++i)
	{
		size_t level = tree[i].get_level();
		location_t loc = tree[i].get_bounding_box().get_center();
	}
	//! [vector traversing]

	//! [level traversing]
	for (size_t i = tree.level_begin(2); i < tree.level_end(2); ++i)
		cluster_t const &c = tree[i];
	//! [level traversing]

	//! [leaf traversing]
	for (auto i : tree.get_leaf_src_indices())
		cluster_t const &c = tree[i];
	//! [leaf traversing]

	//! [source receiver traversing]
	for (auto i : tree.get_leaf_src_indices())
		cluster_t const &c = tree[i];
	//! [source receiver traversing]

	delete[] begin;
	delete[] src_begin;
	delete[] rec_begin;

	return 0;
}
