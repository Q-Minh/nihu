#ifndef LSET_HPP_INCLUDED
#define LSET_HPP_INCLUDED

#include "../tmp/integer.hpp"
#include "domain.hpp"

template <class lset>
struct lset_traits;

template <class Derived>
class LSetBase
{
public:
	typedef lset_traits<Derived> traits;
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



struct line_1_lset;		// linear line lset (-1; +1)
struct tria_1_lset;		// linear triangle lset
struct quad_1_lset;		// linear quad lset
struct brick_1_lset;	// linear brick lset

template<>
struct lset_traits<line_1_lset>
{
	typedef line_domain	domain;
	typedef int_<2>		num_nodes;
};

template<>
struct lset_traits<tria_1_lset>
{
	typedef tria_domain	domain;
	typedef int_<3>		num_nodes;
};

template<>
struct lset_traits<quad_1_lset>
{
	typedef quad_domain	domain;
	typedef int_<4>		num_nodes;
};

template<>
struct lset_traits<brick_1_lset>
{
	typedef brick_domain	domain;
	typedef int_<8>			num_nodes;
};


class line_1_lset : public LSetBase<line_1_lset>
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


class tria_1_lset : public LSetBase<tria_1_lset>
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
		dL << -1.0, -1.0,
			1.0, 0.0,
			0.0, 1.0;
		return dL;
	}
};


class quad_1_lset : public LSetBase<quad_1_lset>
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


template <class lset_from, class lset_to>
struct lset_converter;

template <>
struct lset_converter<quad_1_lset, tria_1_lset>
{
	template <int dim>
	static bool eval(Matrix<double, 4, dim> const &coords)
	{
		tria_1_lset::xi_type xi;
		xi << 1.0, 1.0;
		tria_1_lset::L_type L = tria_1_lset::eval_L(xi);
		return (L.transpose() * coords.topRows(3) - coords.row(3)).norm() < 1e-3;
	}
};

#endif

