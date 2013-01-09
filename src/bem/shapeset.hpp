#ifndef SHAPESET_HPP_INCLUDED
#define SHAPESET_HPP_INCLUDED

#include "../tmp/integer.hpp"
#include "domain.hpp"

template <class shape_set>
struct shape_set_traits;

template <class Derived>
class ShapeSetBase
{
public:
	typedef shape_set_traits<Derived> traits;
	typedef typename traits::domain domain;
	static int const num_nodes = traits::num_nodes::value;

	typedef typename domain::xi_type xi_type;
	typedef Matrix<double, num_nodes, 1> L_type;
	typedef Matrix<double, num_nodes, domain::dimension> dL_type;

	static L_type eval_L(xi_type const &xi)
	{
		return Derived::eval_L(xi);
	}

	static dL_type eval_dL(xi_type const &xi)
	{
		return Derived::eval_dL(xi);
	}
};


struct line_1_shape_set;		// linear line shape_set (-1; +1)
struct tria_1_shape_set;		// linear triangle shape_set
struct parallelogram_shape_set;		// linear quad shape_set
struct quad_1_shape_set;		// linear quad shape_set
struct brick_1_shape_set;	// linear brick shape_set

template<>
struct shape_set_traits<line_1_shape_set>
{
	typedef line_domain	domain;
	typedef int_<2>		num_nodes;
};

template<>
struct shape_set_traits<tria_1_shape_set>
{
	typedef tria_domain	domain;
	typedef int_<3>		num_nodes;
};

template<>
struct shape_set_traits<parallelogram_shape_set>
{
	typedef quad_domain	domain;
	typedef int_<3>		num_nodes;
};

template<>
struct shape_set_traits<quad_1_shape_set>
{
	typedef quad_domain	domain;
	typedef int_<4>		num_nodes;
};

template<>
struct shape_set_traits<brick_1_shape_set>
{
	typedef brick_domain	domain;
	typedef int_<8>			num_nodes;
};


class line_1_shape_set : public ShapeSetBase<line_1_shape_set>
{
public:
	static L_type eval_L(xi_type const &xi)
	{
		L_type L;
		L <<
			(1.0-xi[0])/2.0,
			(1.0+xi[0])/2.0;
		return L;
	}

	static dL_type eval_dL(xi_type const &)
	{
		dL_type dL;
		dL <<
			-.5,
			+.5;
		return dL;
	}
};


class tria_1_shape_set : public ShapeSetBase<tria_1_shape_set>
{
public:
	static L_type eval_L(xi_type const &xi)
	{
		L_type L;
		L <<
			1.0-xi[0]-xi[1],
			xi[0],
			xi[1];
		return L;
	}

	static dL_type eval_dL(xi_type const &)
	{
		dL_type dL;
		dL <<
			-1.0, -1.0,
			+1.0,  0.0,
			 0.0, +1.0;
		return dL;
	}
};


class parallelogram_shape_set : public ShapeSetBase<parallelogram_shape_set>
{
public:
	static L_type eval_L(xi_type const &xi)
	{
		L_type L;
		L <<
			(-xi[0]-xi[1])/2.0,
			(1.0+xi[0])/2.0,
			(1.0+xi[1])/2.0;
		return L;
	}

	static dL_type eval_dL(xi_type const &xi)
	{
		dL_type dL;
		dL <<
			-1.0, -1.0,
			+1.0,  0.0,
			 0.0, +1.0;
		return dL;
	}
};


class quad_1_shape_set : public ShapeSetBase<quad_1_shape_set>
{
public:
	static L_type eval_L(xi_type const &xi)
	{
		L_type L;
		L <<
			(1.0-xi[0])*(1.0-xi[1])/4.0,
			(1.0+xi[0])*(1.0-xi[1])/4.0,
			(1.0+xi[0])*(1.0+xi[1])/4.0,
			(1.0-xi[0])*(1.0+xi[1])/4.0;
		return L;
	}

	static dL_type eval_dL(xi_type const &xi)
	{
		dL_type dL;
		dL <<
			-.25 * (1-xi[1]), (1.0-xi[0]) * -.25,
			+.25 * (1-xi[1]), (1.0+xi[0]) * -.25,
			+.25 * (1+xi[1]), (1.0+xi[0]) * +.25,
			-.25 * (1+xi[1]), (1.0-xi[0]) * +.25;
		return dL;
	}
};


template <class shape_set_from, class shape_set_to>
struct shape_set_converter;

template <>
struct shape_set_converter<quad_1_shape_set, parallelogram_shape_set>
{
	typedef quad_1_shape_set from_set;
	typedef parallelogram_shape_set to_set;
	
	template <int dim>
	static bool eval(Matrix<double, 4, dim> const &coords)
	{
		Matrix<double, 3, dim> c;
		c.row(0) = coords.row(0);
		c.row(1) = coords.row(1);
		c.row(2) = coords.row(3);
		
		to_set::xi_type xi;
		xi << 1.0, 1.0;
		
		return (to_set::eval_L(xi).transpose() * c - coords.row(2)).norm() < 1e-3;
	}
};

#endif

