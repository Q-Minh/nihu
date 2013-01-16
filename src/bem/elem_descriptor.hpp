#ifndef ELEMDESCRIPTOR_HPP_INCLUDED
#define ELEMDESCRIPTOR_HPP_INCLUDED

#include "element.hpp"
#include "quadrature.hpp"
#include <array>

template <class Elem>
class ElemDescriptor
{
public:
	typedef Elem elem_type;
	typedef typename elem_type::x_type x_type;

	ElemDescriptor() {} // default constructor needed for array container in ElemAccelerator
	ElemDescriptor(x_type const &x, x_type const &normal, double jac) : x(x), normal(normal/normal.norm()), jacobian(jacobian) {}

	double get_jacobian(void) const
	{
		return jacobian;
	}

	x_type const &get_x(void) const
	{
		return x;
	}

	x_type const &get_normal(void) const
	{
		return normal;
	}
protected:
	x_type x;
	x_type normal;
	double jacobian;
};

template <class Elem, unsigned N>
class ElemAccelerator
{
public:
	typedef Elem elem_t;
	typedef GaussQuad<typename elem_t::lset::domain, N> quadrature_t;
	static const unsigned size = quadrature_t::size;
	typedef std::array<ElemDescriptor<elem_t>, size> container_t;
	typedef typename container_t::const_iterator iterator_t;

	ElemAccelerator(elem_t const &e)
		: e(e)
	{
		quadrature_t::xi_type const &xi_vec = quadrature_t::get_xi();
		typedef quadrature_t::xi_type::Index index;
		for (index i = 0; i < size; ++i)
		{
			elem_t::x_type x = e.get_x(xi_vec.row(i));
			elem_t::x_type normal = e.get_normal(xi_vec.row(i));
			double jacobian = normal.norm();
			normal /= jacobian;
			data[i] = ElemDescriptor<elem_t>(x, normal, jacobian);
		}
	}

	iterator_t begin(void) const { return data.begin(); }
	iterator_t end(void) const { return data.end(); }

protected:
	elem_t const &e;
	container_t data;
};

#endif
