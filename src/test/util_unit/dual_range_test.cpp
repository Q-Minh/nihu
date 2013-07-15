#include "util/dual_range.hpp"
#include <iostream>

template <class range>
void iterate(range r)
{
	for (auto it = r.begin(); it != r.end(); ++it)
		std::cout << *(it.get_first()) << ' ' << *(it.get_second()) << std::endl;
}

int main(void)
{
	int a[] = {0, 1, 2, 3};
	char b[] = {'a', 'b', 'c', 'd'};
	
	std::cout << "diagonal iteration" << std::endl;
	iterate(create_dual_range(iteration::diagonal(), a, a+4, b, b+4));

	std::cout << "matrix iteration" << std::endl;
	iterate(create_dual_range(iteration::plain(), a, a+4, b, b+4));

	std::cout << "matrix iteration with empty inner" << std::endl;
	iterate(create_dual_range(iteration::plain(), a, a+4, b, b));

	std::cout << "matrix iteration with empty outer" << std::endl;
	iterate(create_dual_range(iteration::plain(), a, a, b, b+4));
	
	return 0;
}

