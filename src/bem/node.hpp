#ifndef NODE_HPP
#define NODE_HPP

#include <Eigen/Dense>
using Eigen::Matrix;

template <int nDim>
class Coord : public Matrix<double, nDim, 1>
{
	typedef Coord type; /* self-returning structure */
};

#endif

