#include "quadrature.hpp"
#include <iostream>

int main(void)
{
	typedef gauss_quad<line_domain, 10> G;
	G::init();

	std::cout << G::get_xi() << std::endl << std::endl << G::get_weight() << std::endl;
	
	return 0;
}

