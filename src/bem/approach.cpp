/*
#include "../tmp/vector.hpp"
#include "kernel.hpp"
#include "function_space.hpp"
*/

#include "vararg.hpp"
#include <iostream>

/*
template <class A> class dirac;

template <class A, class B> class bind_space;

typedef green_HG_kernel G_t;
typedef tmp::vector<tria_1_elem, quad_1_elem> elements_t;
typedef Mesh<elements_t> mesh_t;
typedef function_space<mesh_t, isoparametric_field> Ny_t;
typedef dirac<Ny_t> Nx_t;
typedef bind_space<bind_space<G_t, Ny_t>, Nx_t> bem_t;
*/

using namespace std;

int main(void)
{
	vararg<int, double, char> b(3, 2.5, 'a');
	std::cout << b.get<0>() << b.get<1>() << b.get<2>() << std::endl;

	return 0;
}
