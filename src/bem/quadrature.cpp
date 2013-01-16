#include "quadrature.hpp"

#include <iostream>


int main(void)
{
	typedef GaussQuad<line_domain, 10> G;
	G::init();

	std::cout << G::get_xi() << std::endl << std::endl << G::get_weight() << std::endl;
	
	return 0;
}

