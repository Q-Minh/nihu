#ifndef DOMAIN_HPP_INCLUDED
#define DOMAIN_HPP_INCLUDED

#include "../tmp/integer.hpp"

#include <Eigen/Dense>
using Eigen::Matrix;

struct line_domain;
struct tria_domain;
struct quad_domain;
struct brick_domain;

template <class domain>
struct domain_traits;

template <>
struct domain_traits<line_domain>
{
	typedef int_<1>	dimension;
};

template <>
struct domain_traits<tria_domain>
{
	typedef int_<2>	dimension;
};

template <>
struct domain_traits<quad_domain>
{
	typedef int_<2>	dimension;
};

template <>
struct domain_traits<brick_domain>
{
	typedef int_<3>	dimension;
};

template <class Derived>
class DomainBase 
{
public:
	typedef domain_traits<Derived> traits;
	static int const dimension = traits::dimension::value;

	typedef Matrix<double, dimension, 1> xi_type;
};

struct line_domain : public DomainBase<line_domain>
{
};

struct tria_domain : public DomainBase<tria_domain>
{
};

struct quad_domain : public DomainBase<quad_domain>
{
};

struct brick_domain : public DomainBase<brick_domain>
{
};

#endif

