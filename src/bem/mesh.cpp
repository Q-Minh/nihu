#include "mesh.hpp"

#include <iostream>
#include <algorithm>


typedef tiny<TriaElem, LinQuadElem, PlaneQuadElem, QuadElem> ElemTypes;

int main(void)
{
	Mesh<3, ElemTypes> mesh;
	
	double node[][3] = {
		{0, 0, 0},
		{1, 0, 0},
		{1, 1, 0},
		{0, 1, 0},
		{0, 0, 1},
		{1, 0, 1},
		{1, 1, 1},
		{0, 1, 2}
	};
	for (size_t i = 0; i < sizeof(node)/sizeof(node[0]); ++i)
		mesh.add_node(node[i]);

	int elem[][4+1] = {
		{3, 0, 1, 2, 0},
		{3, 0, 2, 3, 0},
		{4, 4, 5, 6, 7},
		{4, 0, 1, 5, 4},
		{4, 1, 2, 6, 5},
		{4, 2, 3, 7, 6},
		{4, 3, 0, 4, 7}
	};
	for (size_t i = 0; i < sizeof(elem)/sizeof(elem[0]); ++i)
		mesh.add_elem(elem[i]);
		
	std::cout << "Tria elements" << std::endl;
	std::for_each(
		mesh.elements.std::vector<TriaElem>::begin(),
		mesh.elements.std::vector<TriaElem>::end(),
		[] (TriaElem &e) { std::cout << e.transpose() << std::endl; }
	);
		
	std::cout << "LinQuad elements" << std::endl;
	std::for_each(
		mesh.elements.std::vector<LinQuadElem>::begin(),
		mesh.elements.std::vector<LinQuadElem>::end(),
		[] (LinQuadElem &e) { std::cout << e.transpose() << std::endl; }
	);
		
	std::cout << "PlaneQuad elements" << std::endl;
	std::for_each(
		mesh.elements.std::vector<PlaneQuadElem>::begin(),
		mesh.elements.std::vector<PlaneQuadElem>::end(),
		[] (PlaneQuadElem &e) { std::cout << e.transpose() << std::endl; }
	);
		
	std::cout << "Quad elements" << std::endl;
	std::for_each(
		mesh.elements.std::vector<QuadElem>::begin(),
		mesh.elements.std::vector<QuadElem>::end(),
		[] (QuadElem &e) { std::cout << e.transpose() << std::endl; }
	);
		

	return 0;
}


