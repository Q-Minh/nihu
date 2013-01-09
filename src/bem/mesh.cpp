#include "mesh.hpp"

typedef tiny<TriaElem, LinQuadElem, QuadElem> ElemTypes;

int main(void)
{
	unsigned e[] = {
		3, 1, 2, 3, 0, 0,
		4, 4, 5, 6, 7, 0,
		3, 1, 2, 3, 0, 0,
		5, 1, 2, 3, 4, 5,
		4, 4, 5, 6, 7, 0
	};

	Mesh<3, ElemTypes> mesh;

	mesh.add_elem(e);
	mesh.add_elem(e+6);
	mesh.add_elem(e+12);
	mesh.add_elem(e+18);
	mesh.add_elem(e+24);

	return 0;
}

