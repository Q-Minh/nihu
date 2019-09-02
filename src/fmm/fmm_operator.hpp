/**
 * @file fmm_operator.hpp 
 * @brief FMM operator types and tags 
 * @ingroup fmm_ops
 */

#ifndef FMM_OPERATOR_HPP_INCLUDED
#define FMM_OPERATOR_HPP_INCLUDED

namespace NiHu
{
namespace fmm
{

struct l2l_tag {};
struct m2m_tag {};
struct m2l_tag {};

struct m2p_tag {};
struct l2p_tag {};
struct p2l_tag {};
struct p2m_tag {};

struct p2p_tag {};

/** @brief Operator defining its tag type */
template <class FmmTag>
struct fmm_operator
{
	typedef FmmTag fmm_tag;
};

} // end of namespace fmm
} // end of namespace NiHu


#endif
