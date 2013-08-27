#include "tmp/stack.hpp"
#include "tmp/algorithm.hpp"

#include <iostream>

int main(void)
{
	typedef tmp::empty_stack myStack_0;
	typedef tmp::push_front<myStack_0, tmp::int_<1> >::type myStack_1;
	typedef tmp::push_front<myStack_1, tmp::int_<2> >::type myStack_2;
	typedef tmp::pop_front<myStack_2>::type myStack_3;
	
	std::cout << "Testing stack and operations" << std::endl;
	std::cout << "============================" << std::endl;
	std::cout << "size(<>)                  = " << tmp::size<myStack_0>::type::value << std::endl;
	std::cout << "size(push(<>))            = " << tmp::size<myStack_1>::type::value << std::endl;
	std::cout << "size(push(push(<>)))      = " << tmp::size<myStack_2>::type::value << std::endl;
	std::cout << "size(pop(push(push(<>)))) = " << tmp::size<myStack_3>::type::value << std::endl;
	std::cout << std::endl;
	
	return 0;
}

