#ifndef ELEMDESCRIPTOR_HPP_INCLUDED
#define ELEMDESCRIPTOR_HPP_INCLUDED

#include "element.hpp"

template <class xType>
class Descriptor
{
public:
	typedef xType x_t;

	/**
	 * \brief default constructor needed for array container in ElemAccelerator
	 */
	Descriptor()
	{
	}
	
	Descriptor(x_t const &x, x_t const &normal, double jacobian)
		: x(x), normal(normal/normal.norm()), jacobian(jacobian)
	{
	}

	double get_jacobian(void) const
	{
		return jacobian;
	}

	x_t const &get_x(void) const
	{
		return x;
	}

	x_t const &get_normal(void) const
	{
		return normal;
	}
	
protected:
	x_t x;
	x_t normal;
	double jacobian;
};


#include "quadrature.hpp"
#include <array>

template <class Elem, unsigned N>
class ElemAccelerator
{
public:
	typedef Elem elem_t;
	typedef typename elem_t::x_t x_t;
	typedef typename elem_t::lset_t::domain_t domain_t;
	typedef gauss_quad<domain_t, N> quadrature_t;
	static const unsigned size = quadrature_t::size;
	typedef std::array<Descriptor<x_t>, size> container_t;
	typedef typename container_t::const_iterator iterator_t;
	typedef typename quadrature_t::xivec_t xivec_t;
	typedef typename quadrature_t::weightvec_t weightvec_t;

	ElemAccelerator(elem_t const &e) : e(e)
	{
		xivec_t const &xi_vec = quadrature_t::get_xi();
		weightvec_t const &weight_vec = quadrature_t::get_weight();
		typedef typename xivec_t::Index index;
		for (index i = 0; i < size; ++i)
		{
			x_t x = e.get_x(xi_vec.row(i));
			x_t normal = e.get_normal(xi_vec.row(i));
			double jacobian = normal.norm();
			normal /= jacobian;
			data[i] = Descriptor<x_t>(x, normal, jacobian*weight_vec[i]);
		}
	}

	iterator_t begin(void) const
	{
		return data.begin();
	}
	
	iterator_t end(void) const
	{
		return data.end();
	}

protected:
	elem_t const &e;
	container_t data;
};

#endif

