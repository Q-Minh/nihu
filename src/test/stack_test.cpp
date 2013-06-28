
#include "stack.hpp"
#include "algorithm.hpp"
using namespace tmp;

#include <iostream>

int main(void)
{
	typedef empty_stack myStack_0;
	typedef push_front<myStack_0, int_<1> >::type myStack_1;
	typedef push_front<myStack_1, int_<2> >::type myStack_2;
	
	std::cout << size<myStack_2>::type::value;
	
	return 0;
}

