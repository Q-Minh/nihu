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

template <class FmmTag>
struct fmm_operator
{
	typedef FmmTag fmm_tag;
};

}

}


#endif
