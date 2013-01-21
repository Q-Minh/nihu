#ifndef DESCRIPTOR_HPP_INCLUDED
#define DESCRIPTOR_HPP_INCLUDED

#include "element.hpp"
#include "../tmp/general.hpp"

template <class xType>
class location
{
public:
	typedef xType x_t;

	template <class elem_t>
	location(elem_t const &elem, typename elem_t::xi_t const &xi)
	{
		static_assert(tmp::is_same<x_t, typename elem_t::x_t>::value, "Element and descriptor location types must match");
		x = elem.get_x(xi);
		jacobian = elem.get_normal(xi).norm();
	}

	double get_jacobian(void) const
	{
		return jacobian;
	}

	x_t const &get_x(void) const
	{
		return x;
	}

protected:
	x_t x;
	double jacobian;
};


#endif

