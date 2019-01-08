#include "library/lib_element.hpp"

int main(void)
{
	NiHu::tria_1_elem::xi_t xi = NiHu::tria_1_elem::xi_t::Random();
	auto L1 = NiHu::tria_1_shape_set::eval_shape<0>(xi);
	auto L2 = NiHu::tria_1_shape_set::eval_shape<0>(xi.reverse());
	
	NiHu::tria_1_elem::coords_t coords = NiHu::tria_1_elem::coords_t::Random();
	NiHu::tria_1_elem::nodes_t nodes = NiHu::tria_1_elem::nodes_t::Random();
	NiHu::tria_1_elem e1(coords, 0, nodes);
	NiHu::tria_1_elem e2(coords.reverse(), 0, nodes.reverse());
	
	auto f2 = e2.flip();
	
	(void)L1;
	(void)L2;
	(void)f2;
	
	return 0;
}
