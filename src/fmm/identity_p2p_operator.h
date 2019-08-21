/// \file identity_p2p_operator.h
/// \definition of class NiHu::fmm::identity_p2p_operator

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

}
}

#endif
