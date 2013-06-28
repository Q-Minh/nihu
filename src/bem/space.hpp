#ifndef SPACE_HPP_INCLUDED
#define SPACE_HPP_INCLUDED

#include "includes.h"

template <class Scalar, unsigned Dimension>
class space
{
public:
	typedef Scalar scalar_t;
	static unsigned const dimension = Dimension;

	typedef Eigen::Matrix<scalar_t, 1, Dimension> location_t;
};

typedef space<double, 2> space_2d;
typedef space<double, 3> space_3d;

#endif // SPACE_HPP_INCLUDED
