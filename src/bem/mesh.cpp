#include "mesh.hpp"

typedef tiny<TriaElem, QuadElem> ElemTypes;

int main(void)
{
	unsigned e[] = {
		3, 1, 2, 3, 0,
		4, 4, 5, 6, 7
	};

	Mesh<3, ElemTypes> mesh;

	mesh.add_elem(e);
	mesh.add_elem(e+5);

	return 0;
}

