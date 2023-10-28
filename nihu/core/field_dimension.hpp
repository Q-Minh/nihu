#ifndef FIELD_DIMENSION_HPP_INCLUDED
#define FIELD_DIMENSION_HPP_INCLUDED

#include <type_traits>

namespace NiHu
{

namespace field_dimension
{

/** \brief type indicating a 1D valued function space */
typedef std::integral_constant<unsigned, 1> _1d;
/** \brief type indicating a 2D valued function space */
typedef std::integral_constant<unsigned, 2> _2d;
/** \brief type indicating a 3D valued function space */
typedef std::integral_constant<unsigned, 3> _3d;

}

}

#endif /* FIELD_DIMENSION_HPP_INCLUDED */
