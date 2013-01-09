#ifndef ELEMENT_HPP_INCLUDED
#define ELEMENT_HPP_INCLUDED

#include "shapeset.hpp"

template <class element>
struct element_traits;

template <class Derived>
class ElementBase
{
public:
	typedef element_traits<Derived> traits;

	static int const dimension = traits::dimension::value;
	typedef typename traits::lset lset;

	static int const num_nodes = lset::num_nodes;
	typedef typename lset::xi_type xi_type;
	typedef typename lset::L_type L_type;
	typedef typename lset::dL_type dL_type;

	typedef Matrix<double, 1, dimension> x_type;
	typedef Matrix<double, num_nodes, dimension> coords_type;

protected:
	coords_type coords;

public:
	ElementBase(coords_type const &coords) : coords(coords) {}

	x_type get_x(xi_type const &xi)
	{
		return lset::eval_L(xi).transpose() * coords;
	}

};

class tria_1_elem;

template <>
struct element_traits<tria_1_elem>
{
	typedef tria_1_shape_set lset;
	typedef int_<3> dimension;
};

class tria_1_elem : public ElementBase<tria_1_elem>
{
public:
	tria_1_elem(coords_type const &coords) : ElementBase(coords) {}
};

class parallelogram_elem;

template <>
struct element_traits<parallelogram_elem>
{
	typedef parallelogram_shape_set lset;
	typedef int_<3> dimension;
};

class parallelogram_elem : public ElementBase<parallelogram_elem>
{
public:
	parallelogram_elem(coords_type const &coords) : ElementBase(coords) {}
};

class quad_1_elem;

template <>
struct element_traits<quad_1_elem>
{
	typedef quad_1_shape_set lset;
	typedef int_<3> dimension;
};

class quad_1_elem : public ElementBase<quad_1_elem>
{
public:
	quad_1_elem(coords_type const &coords) : ElementBase(coords) {}
};

#endif

