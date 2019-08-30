/** 
 * \file identity_p2p_operator.h
 * \brief Definition of class \ref NiHu::fmm::identity_p2p_operator
 * \ingroup fmm_ops
 */

#ifndef IDENTITY_P2P_OPERATOR_H_INCLUDED
#define IDENTITY_P2P_OPERATOR_H_INCLUDED

#include "fmm_operator.hpp"

namespace NiHu
{
namespace fmm
{

class identity_p2p_operator
	: public fmm_operator<p2p_tag>
{
};

} // end of namespace fmm
} // end of namespace NiHu

#endif /* IDENTITY_P2P_OPERATOR_H_INCLUDED */
