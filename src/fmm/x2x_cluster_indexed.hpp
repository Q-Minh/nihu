/// \file x2x_cluster_indexed.hpp
/// \brief implementation of class fmm::x2x_cluster_indexed
#ifndef FMM_X2X_CLUSTER_INDEXED_HPP_INCLUDED
#define FMM_X2X_CLUSTER_INDEXED_HPP_INCLUDED

#include "cluster_tree.hpp"
#include "fmm_operator.hpp"

#include <type_traits>

namespace NiHu
{
namespace fmm
{

/// \brief class performing indexing of an X2X operator with a cluster
/// \tparam Operator the X2X operator
template <class Operator>
class x2x_cluster_indexed
	: public fmm_operator<typename std::decay<Operator>::type::fmm_tag>
{
public:
	typedef typename std::decay<Operator>::type operator_t;
	typedef typename operator_t::cluster_t cluster_t;
	typedef cluster_tree<cluster_t> tree_t;
	typedef typename operator_t::result_t result_t;

	x2x_cluster_indexed(Operator &&op, tree_t const &tree)
		: m_op(std::forward<Operator>(op))
		, m_tree(tree)
	{
	}

	result_t operator()(size_t to, size_t from) const
	{
		return m_op(m_tree[to], m_tree[from]);
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
x2x_cluster_indexed<Operator>
create_x2x_cluster_indexed(Operator &&op,
	cluster_tree<typename std::decay<Operator>::type::cluster_t> const &tree)
{
	return 	x2x_cluster_indexed<Operator>(std::forward<Operator>(op), tree);
}

} // end of namespace fmm
} // namespace NiHu

#endif //  FMM_X2X_CLUSTER_INDEXED_HPP_INCLUDED
