#include "core/element.hpp"
#include "library/lib_shape.hpp"
#include "core/gaussian_quadrature.hpp"

#include <iostream>

void Quad()
{
//! [define]
	NiHu::gaussian_quadrature<NiHu::line_domain> quadrature(7);
//! [define]

//! [traverse]
	for (auto q : quadrature)
		std::cout << q.get_xi() << '\t' << q.get_w() << '\n';
//! [traverse]

//! [traverse 2]
	for (auto it = quadrature.begin(); it != quadrature.end(); ++it)
		std::cout << it->get_xi() << '\t' << it->get_w() << '\n';
//! [traverse 2]

//! [3D quad]
	NiHu::gaussian_quadrature<NiHu::quad_domain> quad2(7);
	for (auto q : quad2)
		std::cout << q.get_xi().transpose() << '\t' << q.get_w() << '\n';
//! [3D quad]
}

void ElemInt()
{
	typedef NiHu::volume_element<NiHu::quad_1_shape_set, double> Element;
	typedef Element::domain_t Domain;
	typedef NiHu::gaussian_quadrature<Domain> Quadrature;

	Element::coords_t coords;
	coords <<
		0., 1., 1., 0.,
		0., 0., 1., 1.;

	Element elem(coords);
	Quadrature quadrature(5);

	for (auto q : quadrature)
		std::cout << elem.get_x(q.get_xi()).transpose() << std::endl;
}

int main(void)
{
	Quad();
	ElemInt();
	return 0;
}

