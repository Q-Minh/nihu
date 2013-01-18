#ifndef DOMAIN_HPP_INCLUDED
#define DOMAIN_HPP_INCLUDED

#include <Eigen/Dense>
using Eigen::Matrix;

template <unsigned Dimension, unsigned ID>
class Domain
{
public:
	static unsigned const dimension = Dimension;
	static unsigned const id = ID;
	typedef Matrix<double, dimension, 1> xi_t;
};

class line_domain : public Domain<1, 2> {};
class tria_domain : public Domain<2, 3> {};
class quad_domain : public Domain<2, 4> {};
class brick_domain : public Domain<3, 8> {};

#endif

