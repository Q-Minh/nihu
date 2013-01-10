#ifndef DOMAIN_HPP_INCLUDED
#define DOMAIN_HPP_INCLUDED

#include <Eigen/Dense>
using Eigen::Matrix;

template <unsigned Dimension>
class domain
{
public:
	static unsigned const dimension = Dimension;
	typedef Matrix<double, dimension, 1> xi_type;
};

class line_domain : public domain<1> {};
class tria_domain : public domain<2> {};
class quad_domain : public domain<2> {};
class brick_domain : public domain<3> {};

#endif

