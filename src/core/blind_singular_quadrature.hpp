#ifndef BLIND_SINGULAR_QUADRATURE_HPP_INCLUDED
#define BLIND_SINGULAR_QUADRATURE_HPP_INCLUDED

#include "formalism.hpp"
#include "blind_transform_selector.hpp"
#include "duffy_quadrature.hpp"

template <class BlindTransform, class RegularFamily, class LSet>
struct blind_singular_quadrature;

template <class RegularFamily, class LSet>
struct blind_singular_quadrature<blind_transform::duffy, RegularFamily, LSet>
{
	typedef duffy_quadrature<RegularFamily, LSet> type;
};

template <class RegularFamily, class LSet>
struct blind_singular_quadrature<blind_transform::square, RegularFamily, LSet>
{
	typedef duffy_quadrature<RegularFamily, LSet> type;
};


#endif // BLIND_SINGULAR_QUADRATURE_HPP_INCLUDED
