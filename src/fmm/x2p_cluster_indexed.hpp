#ifndef FMM_X2P_CLUSTER_INDEXED_HPP_INCLUDED
#define FMM_X2P_CLUSTER_INDEXED_HPP_INCLUDED

#include "cluster_tree.hpp"
#include "util/matrix_traits.hpp"

#include <type_traits>

namespace NiHu
{
namespace fmm
{

template <class Operator>
class x2p_cluster_indexed
{
public:
	typedef typename std::decay<Operator>::type operator_t;
	typedef typename operator_t::trial_input_t cluster_t;
	typedef cluster_tree<cluster_t> tree_t;
	typedef typename operator_t::result_t op_result_t;
	typedef typename scalar<op_result_t>::type scalar_t;
	typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> result_t;
	static size_t const rows = num_rows<op_result_t>::value;
	static size_t const num_dof_per_rec = rows;

	x2p_cluster_indexed(Operator &&op, tree_t const &tree)
		: m_op(std::forward<Operator>(op))
		, m_tree(tree)
	{
	}

	result_t operator()(size_t to, size_t from) const
	{
		cluster_t const &clus_to = m_tree[to];
		cluster_t const &clus_from = m_tree[from];

		/** \todo delegate cols computation to operator */
		bool first = true;
		size_t cols = 0;

		result_t mat;
		for (size_t ii = 0; ii < clus_to.get_n_rec_nodes(); ++ii)
		{
			size_t i = clus_to.get_rec_node_idx()[ii];
			typename operator_t::result_t res = m_op(i, clus_from);
			if (first)
			{
				cols = res.cols();
				mat.resize(clus_to.get_n_rec_nodes() * rows, cols);
				first = false;
			}
			mat.block(rows * ii, 0, rows, cols) = res;
		}

		return mat;
	}

	result_t operator()(size_t idx) const
	{
		return (*this)(idx, idx);
	}

	tree_t const &get_tree() const
	{
		return m_tree;
	}

private:
	Operator m_op;
	tree_t const &m_tree;
};


template <class Operator>
x2p_cluster_indexed<Operator>
create_x2p_cluster_indexed(Operator &&op,
	cluster_tree<typename std::decay<Operator>::type::trial_input_t> const &tree)
{
	return 	x2p_cluster_indexed<Operator>(std::forward<Operator>(op), tree);
}

} // end of namespace fmm
} // namespace NiHu

#endif //  FMM_X2P_CLUSTER_INDEXED_HPP_INCLUDED
