/// \file weighted_input.hpp
/// \brief implementation of metafunction NiHu::weighted_input

#ifndef WEIGHTED_INPUT_HPP_DEFINED
#define WEIGHTED_INPUT_HPP_DEFINED

#include "location_normal.hpp"

namespace NiHu
{

template <class Input, class Elem>
struct weighted_input;

template <class LSet, class Scalar>
struct weighted_input<location_input<typename surface_element<LSet, Scalar>::space_t>,
	surface_element<LSet, Scalar> >
{
	typedef location_normal_jacobian_input<typename surface_element<LSet, Scalar>::space_t> type;
};

template <class LSet, class Scalar>
struct weighted_input<location_input<typename volume_element<LSet, Scalar>::space_t>,
	volume_element<LSet, Scalar> >
{
	typedef location_volume_jacobian_input<typename volume_element<LSet, Scalar>::space_t> type;
};

template <class LSet, class Scalar>
struct weighted_input<location_normal_jacobian_input<typename surface_element<LSet, Scalar>::space_t>,
	surface_element<LSet, Scalar> >
{
	typedef location_normal_jacobian_input<typename surface_element<LSet, Scalar>::space_t> type;
};

template <class LSet, class Scalar>
struct weighted_input<location_volume_jacobian_input<typename volume_element<LSet, Scalar>::space_t>,
	volume_element<LSet, Scalar> >
{
	typedef location_volume_jacobian_input<typename volume_element<LSet, Scalar>::space_t> type;
};

} // end of namespace NiHu

#endif // WEIGTHED_INPUT DEFINED
