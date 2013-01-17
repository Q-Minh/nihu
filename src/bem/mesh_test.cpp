#include "mesh.hpp"

typedef tiny<
	tiny<tria_1_elem>,
	tiny<parallelogram_elem, quad_1_elem>
> elem_list;

/* metafunction returning a vector's first elements' domain */
template <class Seq>
struct firstdomain
{
	typedef typename deref<typename begin<Seq>::type>::type::domain_t type;
};

typedef transform<
	elem_list,
	inserter<tiny<>, push_back<_1, _2> >,
	firstdomain<_1>
>::type domains_t;

#include <iostream>

int main(void)
{
	std::cout << deref<next<begin<domains_t>::type>::type>::type::dimension << std::endl;

	double c[] = {
		0, 0, 0,
		1, 0, 0,
		1, 1, 0,
		0, 1, 0
	};
	unsigned e[] = {
		4,  0, 1, 2, 3,
		3,  0, 1, 2, 0,
		3,  1, 2, 3, 0
	};

	typedef tiny<tria_1_elem, parallelogram_elem, quad_1_elem> ElemTypes;
	Mesh<ElemTypes> mesh;
	for (unsigned i = 0; i < 4; ++i)
		mesh.add_node(c+i*3);
	for (unsigned i = 0; i < 3; ++i)
		mesh.add_elem(e+i*5);

	return 0;
}
