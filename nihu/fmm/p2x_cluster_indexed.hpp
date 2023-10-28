#ifndef FMM_P2X_CLUSTER_INDEXED_HPP_INCLUDED
#define FMM_P2X_CLUSTER_INDEXED_HPP_INCLUDED

#include "cluster_tree.hpp"
#include "fmm_operator.hpp"
#include "nihu/util/matrix_traits.hpp"


#include <type_traits>

namespace NiHu
{
namespace fmm
{

template <class Operator>
class p2x_cluster_indexed
	: public fmm_operator<typename std::decay<Operator>::type::fmm_tag>
{
public:
	typedef typename std::decay<Operator>::type operator_t;
	typedef typename operator_t::test_input_t cluster_t;
	typedef cluster_tree<cluster_t> tree_t;
	typedef typename operator_t::result_t op_result_t;
	typedef typename scalar<op_result_t>::type scalar_t;
	typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> result_t;
	static size_t const cols = num_cols<op_result_t>::value;
	static size_t const num_dof_per_src = cols;

	p2x_cluster_indexed(Operator &&op, tree_t const &tree)
		: m_op(std::forward<Operator>(op))
		, m_tree(tree)
	{
	}

	result_t operator()(size_t to, size_t from) const
	{
		/** \todo delegate rows computation to operator */
		bool first = true;
		size_t rows = 0;
		cluster_t const &clus_to = m_tree[to];
		cluster_t const &clus_from = m_tree[from];

		result_t mat;
		for (size_t jj = 0; jj < clus_from.get_n_src_nodes(); ++jj)
		{
			size_t j = clus_from.get_src_node_idx()[jj];
			typename operator_t::result_t res = m_op(clus_to, j);
			if (first)
			{
				rows = res.rows();
				mat.resize(rows, clus_from.get_n_src_nodes() * cols);
				first = false;
			}
			mat.block(0, cols * jj, rows, cols) = res;
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
p2x_cluster_indexed<Operator>
create_p2x_cluster_indexed(Operator &&op,
	cluster_tree<typename std::decay<Operator>::type::test_input_t> const &tree)
{
	return 	p2x_cluster_indexed<Operator>(std::forward<Operator>(op), tree);
}

} // end of namespace fmm
} // namespace NiHu

#endif //  FMM_P2X_CLUSTER_INDEXED_HPP_INCLUDED
