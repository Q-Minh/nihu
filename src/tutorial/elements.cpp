#include <iostream>
#include "core/element.hpp"
#include "library/lib_shape.hpp"

void volume(void)
{
//! [Elem type]
	typedef volume_element<quad_1_shape_set, double> Element;
//! [Elem type]

//! [Elem coordinates]
	Element::coords_t coords;
	coords <<
		0.0, 1.0, 1.0, 0.0,
		0.0, 0.0, 1.0, 1.0;

	Element elem(coords);
//! [Elem coordinates]

//! [Local coordinate]
	Element::xi_t xi;
	xi << 0.0, 0.0;
//! [Local coordinate]

//! [location]
	std::cout << "location:\n" << elem.get_x(xi) << std::endl;
//! [location]

//! [gradient]
	std::cout << "gradient:\n" << elem.get_dx(xi) << std::endl;
	std::cout << "jacobian:\n" << elem.get_dx(xi).determinant() << std::endl;
//! [gradient]

//! [second derivative]
	std::cout << "2nd derivative:\n" << elem.get_ddx(xi) << std::endl;
//! [second derivative]
}


void surface(void)
{
//! [surface elem type]
	typedef surface_element<quad_1_shape_set, double> Element;

	Element::coords_t coords;
	coords <<
		0.0, 1.0, 1.0, 0.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;

	Element elem(coords);

	Element::xi_t xi;
	xi << 0.0, 0.0;

	std::cout << "location:\n" << elem.get_x(xi) << std::endl;
//! [surface elem type]

//! [surface gradient]
	std::cout << "gradient:\n" << elem.get_dx(xi) << std::endl;
	std::cout << "normal:\n" << elem.get_normal(xi) << std::endl;
	std::cout << "jacobian:\n" << elem.get_normal(xi).norm() << std::endl;
//! [surface gradient]

//! [surface second derivative]
	std::cout << "2nd derivative:\n" << elem.get_ddx(xi) << std::endl;
//! [surface second derivative]
}

