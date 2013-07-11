#include "../tmp/algorithm.hpp"

#include <iostream>

int main(void)
{
	typedef tmp::vector<std::false_type, std::true_type> v;
	typedef tmp::push_back<v, int>::type v2;
	return 0;
}