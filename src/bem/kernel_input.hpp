#ifndef KERNELINPUT_HPP_INCLUDED
#define KERNELINPUT_HPP_INCLUDED

#include "element.hpp"

template <class Element>
class Location
{
public:
	typedef Element elem_type;
	typedef typename elem_type::x_type x_type;
	
	Location(x_type const &x) : x(x) {}
	x_type const &get_x(void) const { return x; }

protected:
	x_type x;
};


template <class Element>
class LocationNormal : public Location<Element>
{
public:
	typedef Element elem_type;
	typedef typename elem_type::x_type x_type;

	LocationNormal(x_type const &x, x_type const &normal) : Location<Element>(x), normal(normal) {}
	x_type const &get_normal(void) const { return normal; }

protected:
	x_type normal;
};

#endif

